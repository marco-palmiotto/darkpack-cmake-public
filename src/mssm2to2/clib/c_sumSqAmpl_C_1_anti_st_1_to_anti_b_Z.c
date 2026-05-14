#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_st_1_to_anti_b_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_st_1_to_anti_b_Z(
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
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_C_1, 2);
    const ccomplex_t IT_0002 = pow(m_st_1, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0004*IT_0006*IT_0008*conj(U_st_00);
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*IT_0007*IT_0012;
    const ccomplex_t IT_0014 = 0.333333333333333*IT_0013;
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = IT_0003*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0007*IT_0012;
    const ccomplex_t IT_0020 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0008*IT_0011;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0021 + -IT_0023 + (-2)
      *IT_0025);
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0010*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0004*IT_0006*IT_0008*conj(U_st_00);
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0035 = IT_0019*IT_0034;
    const ccomplex_t IT_0036 = IT_0022*IT_0034;
    const ccomplex_t IT_0037 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0038 = IT_0022*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + -IT_0036 + (-2)
      *IT_0038);
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = IT_0033*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0031 + IT_0044;
    const ccomplex_t IT_0046 = (-2)*IT_0031 + (-2)*IT_0044;
    const ccomplex_t IT_0047 = s_13*s_14;
    const ccomplex_t IT_0048 = s_34*IT_0001;
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = IT_0047 + IT_0049;
    const ccomplex_t IT_0051 = 6*IT_0050;
    const ccomplex_t IT_0052 = IT_0046*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0008
      *IT_0011 + -1./3*IT_0007*IT_0012)*U_st_00*conj(U_st_00) + IT_0007*IT_0012
      *U_st_10*conj(U_st_10));
    const ccomplex_t IT_0054 = 1.33333333333333*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0010*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0054*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0008
      *IT_0011 + -1./3*IT_0007*IT_0012)*conj(U_st_00)*U_st_01 + IT_0007*IT_0012
      *conj(U_st_10)*U_st_11);
    const ccomplex_t IT_0060 = 1.33333333333333*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0004*IT_0006*IT_0008*conj(U_st_01);
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0060*IT_0065;
    const ccomplex_t IT_0067 = IT_0058 + IT_0066;
    const ccomplex_t IT_0068 = (-6)*IT_0050;
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (-2)*IT_0017 + -IT_0058 + -IT_0066;
    const ccomplex_t IT_0071 = s_14*IT_0000;
    const ccomplex_t IT_0072 = s_13*s_34;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = IT_0071 + IT_0073;
    const ccomplex_t IT_0075 = (-6)*IT_0074;
    const ccomplex_t IT_0076 = IT_0070*IT_0075;
    const ccomplex_t IT_0077 = pow(m_Z, -2);
    const ccomplex_t IT_0078 = pow(m_Z, 4);
    const ccomplex_t IT_0079 = s_13*IT_0078;
    const ccomplex_t IT_0080 = IT_0077*IT_0079;
    const ccomplex_t IT_0081 = 6*IT_0080;
    const ccomplex_t IT_0082 = s_14*s_34;
    const ccomplex_t IT_0083 = 24*IT_0082;
    const ccomplex_t IT_0084 = pow(m_Z, 2);
    const ccomplex_t IT_0085 = s_13*IT_0084;
    const ccomplex_t IT_0086 = (-12)*IT_0085;
    const ccomplex_t IT_0087 = IT_0081 + IT_0083 + IT_0086;
    const ccomplex_t IT_0088 = (-6)*IT_0080;
    const ccomplex_t IT_0089 = IT_0083 + IT_0088;
    const ccomplex_t IT_0090 = conj(IT_0046)*IT_0051;
    const ccomplex_t IT_0091 = sin(beta);
    const ccomplex_t IT_0092 = cpow(IT_0091, -1);
    const ccomplex_t IT_0093 = IT_0008*IT_0092;
    const ccomplex_t IT_0094 = m_t*V_tb*V_u1*e_em*IT_0004*conj(U_st_10);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0100 = IT_0008*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = IT_0098 + IT_0103;
    const ccomplex_t IT_0105 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0106 = IT_0019*IT_0105;
    const ccomplex_t IT_0107 = IT_0022*IT_0105;
    const ccomplex_t IT_0108 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0109 = IT_0022*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0106 + -IT_0107 + (-2)
      *IT_0109);
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0029*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = IT_0104*IT_0113;
    const ccomplex_t IT_0115 = m_t*V_tb*V_u2*e_em*IT_0004*conj(U_st_10);
    const ccomplex_t IT_0116 = IT_0093*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = V_tb*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0121 = IT_0008*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = IT_0119 + IT_0124;
    const ccomplex_t IT_0126 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0127 = IT_0019*IT_0126;
    const ccomplex_t IT_0128 = IT_0022*IT_0126;
    const ccomplex_t IT_0129 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0130 = IT_0022*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + -IT_0128 + (-2)
      *IT_0130);
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = IT_0042*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0125*IT_0134;
    const ccomplex_t IT_0136 = -IT_0114 + -IT_0135;
    const ccomplex_t IT_0137 = m_b*m_C_1;
    const ccomplex_t IT_0138 = IT_0084*IT_0137;
    const ccomplex_t IT_0139 = (-3)*IT_0138;
    const ccomplex_t IT_0140 = 6*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0142 = IT_0027*IT_0104;
    const ccomplex_t IT_0143 = IT_0029*IT_0141*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0145 = IT_0040*IT_0125;
    const ccomplex_t IT_0146 = IT_0042*IT_0144*IT_0145;
    const ccomplex_t IT_0147 = m_C_1*IT_0102;
    const ccomplex_t IT_0148 = m_C_1*IT_0098;
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = IT_0147 + IT_0149;
    const ccomplex_t IT_0151 = IT_0113*IT_0150;
    const ccomplex_t IT_0152 = m_C_1*IT_0123;
    const ccomplex_t IT_0153 = m_C_1*IT_0119;
    const ccomplex_t IT_0154 = -IT_0153;
    const ccomplex_t IT_0155 = IT_0152 + IT_0154;
    const ccomplex_t IT_0156 = IT_0134*IT_0155;
    const ccomplex_t IT_0157 = -IT_0143 + -IT_0146 + IT_0151 + IT_0156;
    const ccomplex_t IT_0158 = e_em*IT_0019;
    const ccomplex_t IT_0159 = e_em*IT_0022;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*(IT_0158 + 3*IT_0159);
    const ccomplex_t IT_0161 = (-0.166666666666667)*IT_0160;
    const ccomplex_t IT_0162 = IT_0010*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0164 = IT_0003*IT_0162*IT_0163;
    const ccomplex_t IT_0165 = m_b*IT_0017;
    const ccomplex_t IT_0166 = IT_0164 + -IT_0165;
    const ccomplex_t IT_0167 = conj(IT_0157) + conj(IT_0166);
    const ccomplex_t IT_0168 = s_14*IT_0084;
    const ccomplex_t IT_0169 = IT_0077*IT_0168;
    const ccomplex_t IT_0170 = 0.5*IT_0169;
    const ccomplex_t IT_0171 = s_14 + IT_0170;
    const ccomplex_t IT_0172 = m_b*IT_0171;
    const ccomplex_t IT_0173 = 2*IT_0172;
    const ccomplex_t IT_0174 = 6*IT_0173;
    const ccomplex_t IT_0175 = IT_0033*IT_0132;
    const ccomplex_t IT_0176 = IT_0042*IT_0144*IT_0175;
    const ccomplex_t IT_0177 = m_C_1*IT_0044;
    const ccomplex_t IT_0178 = IT_0010*IT_0111;
    const ccomplex_t IT_0179 = IT_0029*IT_0141*IT_0178;
    const ccomplex_t IT_0180 = m_C_1*IT_0031;
    const ccomplex_t IT_0181 = IT_0176 + IT_0177 + IT_0179 + IT_0180;
    const ccomplex_t IT_0182 = IT_0014*IT_0104;
    const ccomplex_t IT_0183 = IT_0003*IT_0163*IT_0182;
    const ccomplex_t IT_0184 = IT_0003*IT_0161;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = m_b*IT_0102;
    const ccomplex_t IT_0187 = m_b*IT_0098;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = IT_0186 + IT_0188;
    const ccomplex_t IT_0190 = IT_0185*IT_0189;
    const ccomplex_t IT_0191 = IT_0183 + IT_0190;
    const ccomplex_t IT_0192 = conj(IT_0181) + conj(IT_0191);
    const ccomplex_t IT_0193 = s_34*IT_0084;
    const ccomplex_t IT_0194 = IT_0077*IT_0193;
    const ccomplex_t IT_0195 = (-0.25)*IT_0194;
    const ccomplex_t IT_0196 = s_34 + IT_0195;
    const ccomplex_t IT_0197 = m_C_1*IT_0196;
    const ccomplex_t IT_0198 = (-4)*IT_0197;
    const ccomplex_t IT_0199 = 6*IT_0198;
    const ccomplex_t IT_0200 = conj(IT_0067)*IT_0068;
    const ccomplex_t IT_0201 = conj(IT_0070)*IT_0075;
    const ccomplex_t IT_0202 = (IT_0098 + -IT_0102)*IT_0185;
    const ccomplex_t IT_0203 = -IT_0202;
    const ccomplex_t IT_0204 = 3*IT_0138;
    const ccomplex_t IT_0205 = (-6)*IT_0204;
    const ccomplex_t IT_0206 = pow(s_14, 2);
    const ccomplex_t IT_0207 = IT_0077*IT_0206;
    const ccomplex_t IT_0208 = -IT_0001;
    const ccomplex_t IT_0209 = IT_0207 + IT_0208;
    const ccomplex_t IT_0210 = m_b*IT_0209;
    const ccomplex_t IT_0211 = 6*IT_0210;
    const ccomplex_t IT_0212 = IT_0046*IT_0211;
    const ccomplex_t IT_0213 = s_14*s_34*IT_0077;
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = s_13 + IT_0214;
    const ccomplex_t IT_0216 = m_b*IT_0215;
    const ccomplex_t IT_0217 = 6*IT_0216;
    const ccomplex_t IT_0218 = IT_0070*IT_0217;
    const ccomplex_t IT_0219 = (-6)*IT_0210;
    const ccomplex_t IT_0220 = IT_0067*IT_0219;
    const ccomplex_t IT_0221 = s_13*IT_0209;
    const ccomplex_t IT_0222 = 6*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0224 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0225 = IT_0008*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = IT_0060*IT_0227;
    const ccomplex_t IT_0229 = m_t*V_tb*V_u1*e_em*IT_0004*conj(U_st_11);
    const ccomplex_t IT_0230 = IT_0093*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = (-0.5)*IT_0232;
    const ccomplex_t IT_0234 = IT_0060*IT_0233;
    const ccomplex_t IT_0235 = -IT_0234;
    const ccomplex_t IT_0236 = IT_0228 + IT_0235;
    const ccomplex_t IT_0237 = IT_0223*IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0239 = IT_0054*IT_0102;
    const ccomplex_t IT_0240 = IT_0054*IT_0098;
    const ccomplex_t IT_0241 = -IT_0240;
    const ccomplex_t IT_0242 = IT_0239 + IT_0241;
    const ccomplex_t IT_0243 = IT_0238*IT_0242;
    const ccomplex_t IT_0244 = -IT_0237 + -IT_0243;
    const ccomplex_t IT_0245 = IT_0137*IT_0209;
    const ccomplex_t IT_0246 = 6*IT_0245;
    const ccomplex_t IT_0247 = m_C_1*IT_0215;
    const ccomplex_t IT_0248 = 6*IT_0247;
    const ccomplex_t IT_0249 = (-6)*IT_0221;
    const ccomplex_t IT_0250 = s_13*IT_0215;
    const ccomplex_t IT_0251 = (-6)*IT_0250;
    const ccomplex_t IT_0252 = 2*IT_0114 + 2*IT_0135;
    const ccomplex_t IT_0253 = (-6)*IT_0245;
    const ccomplex_t IT_0254 = IT_0104*IT_0185;
    const ccomplex_t IT_0255 = IT_0237 + IT_0243 + (-2)*IT_0254;
    const ccomplex_t IT_0256 = IT_0137*IT_0215;
    const ccomplex_t IT_0257 = (-6)*IT_0256;
    const ccomplex_t IT_0258 = conj(IT_0136) + conj(IT_0203);
    const ccomplex_t IT_0259 = 6*IT_0250;
    const ccomplex_t IT_0260 = 6*IT_0256;
    const ccomplex_t IT_0261 = (-6)*IT_0247;
    const ccomplex_t IT_0262 = pow(s_34, 2);
    const ccomplex_t IT_0263 = IT_0077*IT_0262;
    const ccomplex_t IT_0264 = -IT_0263;
    const ccomplex_t IT_0265 = IT_0000 + IT_0264;
    const ccomplex_t IT_0266 = m_C_1*IT_0265;
    const ccomplex_t IT_0267 = 6*IT_0266;
    const ccomplex_t IT_0268 = s_13*IT_0265;
    const ccomplex_t IT_0269 = (-6)*IT_0268;
    const ccomplex_t IT_0270 = IT_0137*IT_0265;
    const ccomplex_t IT_0271 = (-6)*IT_0270;
    const ccomplex_t IT_0272 = (-2)*IT_0082;
    const ccomplex_t IT_0273 = IT_0085 + IT_0272;
    const ccomplex_t IT_0274 = IT_0077*IT_0273;
    const ccomplex_t IT_0275 = (-6)*IT_0274;
    const ccomplex_t IT_0276 = 12*s_13;
    const ccomplex_t IT_0277 = IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = conj(IT_0046)*IT_0211;
    const ccomplex_t IT_0279 = conj(IT_0070)*IT_0217;
    const ccomplex_t IT_0280 = conj(IT_0244)*IT_0248;
    const ccomplex_t IT_0281 = conj(IT_0255)*IT_0267;
    const ccomplex_t IT_0282 = (-18)*IT_0137;
    const ccomplex_t IT_0283 = conj(IT_0181)*IT_0282;
    const ccomplex_t IT_0284 = conj(IT_0067)*IT_0219;
    const ccomplex_t IT_0285 = conj(IT_0252)*IT_0261;
    const ccomplex_t IT_0286 = (-0.25)*IT_0169;
    const ccomplex_t IT_0287 = s_14 + IT_0286;
    const ccomplex_t IT_0288 = m_b*IT_0287;
    const ccomplex_t IT_0289 = (-4)*IT_0288;
    const ccomplex_t IT_0290 = (-6)*IT_0289;
    const ccomplex_t IT_0291 = 0.5*IT_0194;
    const ccomplex_t IT_0292 = s_34 + IT_0291;
    const ccomplex_t IT_0293 = m_C_1*IT_0292;
    const ccomplex_t IT_0294 = 2*IT_0293;
    const ccomplex_t IT_0295 = (-6)*IT_0294;
    const ccomplex_t IT_0296 = (IT_0157 + IT_0166)*(conj(IT_0045)*IT_0174 +
       conj(IT_0136)*IT_0199 + IT_0167*IT_0277 + IT_0278 + IT_0279 + IT_0280 +
       IT_0281 + conj(IT_0191)*IT_0282 + IT_0283 + IT_0284 + IT_0285 + conj
      (IT_0018)*IT_0290 + conj(IT_0203)*IT_0295);
    const ccomplex_t IT_0297 = IT_0181*IT_0277;
    const ccomplex_t IT_0298 = IT_0075*conj(IT_0255);
    const ccomplex_t IT_0299 = IT_0068*conj(IT_0244);
    const ccomplex_t IT_0300 = IT_0211*conj(IT_0252);
    const ccomplex_t IT_0301 = IT_0217*conj(IT_0255);
    const ccomplex_t IT_0302 = IT_0219*conj(IT_0244);
    const ccomplex_t IT_0303 = conj(IT_0018)*IT_0295;
    const ccomplex_t IT_0304 = conj(IT_0046)*IT_0261;
    const ccomplex_t IT_0305 = conj(IT_0067)*IT_0248;
    const ccomplex_t IT_0306 = conj(IT_0070)*IT_0267;
    const ccomplex_t IT_0307 = IT_0167*IT_0282;
    const ccomplex_t IT_0308 = conj(IT_0136)*IT_0174;
    const ccomplex_t IT_0309 = conj(IT_0045)*IT_0199;
    const ccomplex_t IT_0310 = (conj(IT_0018) + conj(IT_0045))*(IT_0052 +
       IT_0069 + IT_0076) + IT_0045*(conj(IT_0045)*IT_0087 + conj(IT_0018)
      *IT_0089 + IT_0090 + conj(IT_0136)*IT_0140 + IT_0167*IT_0174 + IT_0192
      *IT_0199 + IT_0200 + IT_0201 + conj(IT_0203)*IT_0205) + IT_0167*(IT_0212 +
       IT_0218 + IT_0220) + IT_0067*(conj(IT_0067)*IT_0222 + conj(IT_0244)
      *IT_0246 + IT_0192*IT_0248 + conj(IT_0046)*IT_0249 + conj(IT_0070)*IT_0251
       + conj(IT_0252)*IT_0253 + conj(IT_0255)*IT_0257) + IT_0244*(IT_0192
      *IT_0219 + IT_0222*conj(IT_0244) + conj(IT_0067)*IT_0246 + IT_0167*IT_0248
       + IT_0249*conj(IT_0252) + conj(IT_0046)*IT_0253 + IT_0251*conj(IT_0255) +
       conj(IT_0070)*IT_0257 + IT_0068*IT_0258) + IT_0252*(IT_0192*IT_0211 +
       conj(IT_0046)*IT_0246 + conj(IT_0244)*IT_0249 + IT_0222*conj(IT_0252) +
       conj(IT_0067)*IT_0253 + IT_0051*IT_0258 + conj(IT_0255)*IT_0259 + conj
      (IT_0070)*IT_0260 + IT_0167*IT_0261) + IT_0046*(conj(IT_0046)*IT_0222 +
       conj(IT_0067)*IT_0249 + IT_0246*conj(IT_0252) + conj(IT_0244)*IT_0253 +
       conj(IT_0070)*IT_0259 + conj(IT_0255)*IT_0260 + IT_0192*IT_0261) +
       IT_0255*(IT_0192*IT_0217 + conj(IT_0244)*IT_0251 + conj(IT_0067)*IT_0257 
      + IT_0075*IT_0258 + conj(IT_0252)*IT_0259 + conj(IT_0046)*IT_0260 +
       IT_0167*IT_0267 + conj(IT_0255)*IT_0269 + conj(IT_0070)*IT_0271) +
       IT_0070*(conj(IT_0067)*IT_0251 + conj(IT_0244)*IT_0257 + conj(IT_0046)
      *IT_0259 + conj(IT_0252)*IT_0260 + IT_0192*IT_0267 + conj(IT_0070)*IT_0269
       + conj(IT_0255)*IT_0271) + IT_0018*(conj(IT_0018)*IT_0087 + conj(IT_0045)
      *IT_0089 + IT_0090 + IT_0200 + IT_0201 + IT_0140*conj(IT_0203) + conj
      (IT_0136)*IT_0205 + IT_0167*IT_0290 + IT_0192*IT_0295) + IT_0296 + IT_0192
      *IT_0297 + IT_0136*(IT_0087*conj(IT_0136) + conj(IT_0045)*IT_0140 +
       IT_0174*IT_0192 + IT_0167*IT_0199 + IT_0089*conj(IT_0203) + conj(IT_0018)
      *IT_0205 + IT_0051*conj(IT_0252) + IT_0298 + IT_0299) + IT_0203*(IT_0089
      *conj(IT_0136) + conj(IT_0018)*IT_0140 + IT_0087*conj(IT_0203) + conj
      (IT_0045)*IT_0205 + IT_0051*conj(IT_0252) + IT_0192*IT_0290 + IT_0167
      *IT_0295 + IT_0298 + IT_0299) + IT_0181*(conj(IT_0203)*IT_0290 + IT_0300 +
       IT_0301 + IT_0302 + IT_0303 + IT_0304 + IT_0305 + IT_0306 + IT_0307 +
       IT_0308 + IT_0309) + IT_0191*(IT_0192*IT_0277 + conj(IT_0203)*IT_0290 +
       IT_0300 + IT_0301 + IT_0302 + IT_0303 + IT_0304 + IT_0305 + IT_0306 +
       IT_0307 + IT_0308 + IT_0309);
    return create_ccomplex_return(IT_0310);
}

