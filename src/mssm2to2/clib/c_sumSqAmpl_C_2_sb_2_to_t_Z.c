#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_2_to_t_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_2_to_t_Z(
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
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_C_2, 2);
    const ccomplex_t IT_0001 = pow(m_sb_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_12 + -IT_0000 + -IT_0001 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = cpow(IT_0005, -1);
    const ccomplex_t IT_0009 = IT_0003*IT_0008;
    const ccomplex_t IT_0010 = e_em*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*(IT_0007 + (-3)*IT_0010);
    const ccomplex_t IT_0012 = (-0.166666666666667)*IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = m_b*conj(U_d2)*V_tb*e_em*IT_0016*U_sb_11;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0021 = IT_0008*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0019 + (-2)*IT_0021);
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0012*IT_0023;
    const ccomplex_t IT_0025 = IT_0002*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = pow(m_Z, -2);
    const ccomplex_t IT_0028 = pow(m_Z, 4);
    const ccomplex_t IT_0029 = s_13*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = 6*IT_0030;
    const ccomplex_t IT_0032 = s_14*s_34;
    const ccomplex_t IT_0033 = 24*IT_0032;
    const ccomplex_t IT_0034 = pow(m_Z, 2);
    const ccomplex_t IT_0035 = s_13*IT_0034;
    const ccomplex_t IT_0036 = (-12)*IT_0035;
    const ccomplex_t IT_0037 = IT_0031 + IT_0033 + IT_0036;
    const ccomplex_t IT_0038 = pow(m_t, 2);
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0001 + IT_0038 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0040 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0041 = IT_0006*IT_0040;
    const ccomplex_t IT_0042 = IT_0009*IT_0040;
    const ccomplex_t IT_0043 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0044 = IT_0009*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0041 + -IT_0042 + (-2)
      *IT_0044);
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = m_b*conj(U_d1)*V_tb*e_em*IT_0016*U_sb_11;
    const ccomplex_t IT_0048 = IT_0015*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0051 = IT_0008*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0049 + (-2)*IT_0051);
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = IT_0046*IT_0053;
    const ccomplex_t IT_0055 = IT_0039*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_23 + IT_0001 + IT_0038 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0058 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0059 = IT_0006*IT_0058;
    const ccomplex_t IT_0060 = IT_0009*IT_0058;
    const ccomplex_t IT_0061 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0062 = IT_0009*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + -IT_0060 + (-2)
      *IT_0062);
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0023*IT_0064;
    const ccomplex_t IT_0066 = IT_0057*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0056 + -IT_0067;
    const ccomplex_t IT_0069 = (-24)*IT_0032;
    const ccomplex_t IT_0070 = IT_0031 + IT_0069;
    const ccomplex_t IT_0071 = s_34*IT_0034;
    const ccomplex_t IT_0072 = IT_0027*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = s_34 + IT_0073;
    const ccomplex_t IT_0075 = m_C_2*IT_0074;
    const ccomplex_t IT_0076 = 2*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0078 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0079 = IT_0006*IT_0078;
    const ccomplex_t IT_0080 = IT_0009*IT_0078;
    const ccomplex_t IT_0081 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0082 = IT_0009*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0079 + -IT_0080 + (-2)
      *IT_0082);
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = IT_0023*IT_0084;
    const ccomplex_t IT_0086 = IT_0057*IT_0077*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0088 = sin(beta);
    const ccomplex_t IT_0089 = cpow(IT_0088, -1);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0008*IT_0016*IT_0089*U_sb_01;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = IT_0012*IT_0091;
    const ccomplex_t IT_0093 = IT_0002*IT_0087*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0095 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0096 = IT_0006*IT_0095;
    const ccomplex_t IT_0097 = IT_0009*IT_0095;
    const ccomplex_t IT_0098 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0099 = IT_0009*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + -IT_0097 + (-2)
      *IT_0099);
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = IT_0053*IT_0101;
    const ccomplex_t IT_0103 = IT_0039*IT_0094*IT_0102;
    const ccomplex_t IT_0104 = m_C_2*IT_0056;
    const ccomplex_t IT_0105 = m_C_2*IT_0067;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0005;
    const ccomplex_t IT_0107 = (-0.666666666666667)*IT_0106;
    const ccomplex_t IT_0108 = IT_0091*IT_0107;
    const ccomplex_t IT_0109 = IT_0002*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = m_t*IT_0110;
    const ccomplex_t IT_0112 = -IT_0086 + -IT_0093 + -IT_0103 + -IT_0104 + 
      -IT_0105 + IT_0111;
    const ccomplex_t IT_0113 = s_14*IT_0034;
    const ccomplex_t IT_0114 = IT_0027*IT_0113;
    const ccomplex_t IT_0115 = (-0.25)*IT_0114;
    const ccomplex_t IT_0116 = s_14 + IT_0115;
    const ccomplex_t IT_0117 = m_t*IT_0116;
    const ccomplex_t IT_0118 = (-4)*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0008*IT_0016*IT_0089*U_sb_01;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = IT_0046*IT_0120;
    const ccomplex_t IT_0122 = IT_0039*IT_0094*IT_0121;
    const ccomplex_t IT_0123 = IT_0064*IT_0091;
    const ccomplex_t IT_0124 = IT_0057*IT_0077*IT_0123;
    const ccomplex_t IT_0125 = IT_0023*IT_0107;
    const ccomplex_t IT_0126 = IT_0002*IT_0087*IT_0125;
    const ccomplex_t IT_0127 = m_t*IT_0026;
    const ccomplex_t IT_0128 = IT_0084*IT_0091;
    const ccomplex_t IT_0129 = IT_0057*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = m_C_2*IT_0130;
    const ccomplex_t IT_0132 = IT_0101*IT_0120;
    const ccomplex_t IT_0133 = IT_0039*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = m_C_2*IT_0134;
    const ccomplex_t IT_0136 = -IT_0122 + -IT_0124 + IT_0126 + -IT_0127 + 
      -IT_0131 + -IT_0135;
    const ccomplex_t IT_0137 = m_t*m_C_2;
    const ccomplex_t IT_0138 = IT_0034*IT_0137;
    const ccomplex_t IT_0139 = (-3)*IT_0138;
    const ccomplex_t IT_0140 = -IT_0110;
    const ccomplex_t IT_0141 = 3*IT_0138;
    const ccomplex_t IT_0142 = -IT_0130 + -IT_0134;
    const ccomplex_t IT_0143 = 6*conj(IT_0026);
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0004
      *IT_0005 + 3*IT_0003*IT_0008)*U_sb_01*conj(U_sb_01) + IT_0004*IT_0005
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0145 = (-0.666666666666667)*IT_0144;
    const ccomplex_t IT_0146 = cpow((-2)*s_13 + IT_0000 + IT_0038 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0147 = IT_0023*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = IT_0145*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0004
      *IT_0005 + 3*IT_0003*IT_0008)*conj(U_sb_00)*U_sb_01 + IT_0004*IT_0005*conj
      (U_sb_10)*U_sb_11);
    const ccomplex_t IT_0151 = (-0.666666666666667)*IT_0150;
    const ccomplex_t IT_0152 = cpow((-2)*s_13 + IT_0000 + IT_0038 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0153 = m_b*conj(U_d2)*V_tb*e_em*IT_0016*U_sb_10;
    const ccomplex_t IT_0154 = IT_0015*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0157 = IT_0008*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0155 + (-2)*IT_0157);
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = IT_0152*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = IT_0151*IT_0161;
    const ccomplex_t IT_0163 = 2*IT_0056 + 2*IT_0067 + -IT_0149 + -IT_0162;
    const ccomplex_t IT_0164 = s_13*s_14;
    const ccomplex_t IT_0165 = s_34*IT_0000;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = IT_0164 + IT_0166;
    const ccomplex_t IT_0168 = 6*IT_0167;
    const ccomplex_t IT_0169 = 2*IT_0026 + IT_0149 + IT_0162;
    const ccomplex_t IT_0170 = s_14*IT_0038;
    const ccomplex_t IT_0171 = s_13*s_34;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = IT_0170 + IT_0172;
    const ccomplex_t IT_0174 = 6*IT_0173;
    const ccomplex_t IT_0175 = 6*IT_0141;
    const ccomplex_t IT_0176 = 6*IT_0139;
    const ccomplex_t IT_0177 = 0.5*IT_0114;
    const ccomplex_t IT_0178 = s_14 + IT_0177;
    const ccomplex_t IT_0179 = m_t*IT_0178;
    const ccomplex_t IT_0180 = 2*IT_0179;
    const ccomplex_t IT_0181 = 6*IT_0180;
    const ccomplex_t IT_0182 = (-0.25)*IT_0072;
    const ccomplex_t IT_0183 = s_34 + IT_0182;
    const ccomplex_t IT_0184 = m_C_2*IT_0183;
    const ccomplex_t IT_0185 = (-4)*IT_0184;
    const ccomplex_t IT_0186 = 6*IT_0185;
    const ccomplex_t IT_0187 = 6*IT_0026;
    const ccomplex_t IT_0188 = IT_0091*IT_0146;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = IT_0145*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0008*IT_0016*IT_0089*U_sb_00;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = IT_0152*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = IT_0151*IT_0194;
    const ccomplex_t IT_0196 = 2*IT_0130 + 2*IT_0134 + IT_0190 + IT_0195;
    const ccomplex_t IT_0197 = (-2)*IT_0110 + -IT_0190 + -IT_0195;
    const ccomplex_t IT_0198 = (-6)*conj(IT_0026);
    const ccomplex_t IT_0199 = 6*IT_0076;
    const ccomplex_t IT_0200 = 6*IT_0118;
    const ccomplex_t IT_0201 = (-6)*IT_0167;
    const ccomplex_t IT_0202 = (-6)*IT_0173;
    const ccomplex_t IT_0203 = pow(s_14, 2);
    const ccomplex_t IT_0204 = IT_0027*IT_0203;
    const ccomplex_t IT_0205 = -IT_0000;
    const ccomplex_t IT_0206 = IT_0204 + IT_0205;
    const ccomplex_t IT_0207 = s_13*IT_0206;
    const ccomplex_t IT_0208 = 6*IT_0207;
    const ccomplex_t IT_0209 = m_t*IT_0206;
    const ccomplex_t IT_0210 = 6*IT_0209;
    const ccomplex_t IT_0211 = IT_0137*IT_0206;
    const ccomplex_t IT_0212 = 6*IT_0211;
    const ccomplex_t IT_0213 = s_14*s_34*IT_0027;
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = s_13 + IT_0214;
    const ccomplex_t IT_0216 = s_13*IT_0215;
    const ccomplex_t IT_0217 = (-6)*IT_0216;
    const ccomplex_t IT_0218 = IT_0137*IT_0215;
    const ccomplex_t IT_0219 = (-6)*IT_0218;
    const ccomplex_t IT_0220 = m_C_2*IT_0215;
    const ccomplex_t IT_0221 = (-6)*IT_0220;
    const ccomplex_t IT_0222 = 12*s_13;
    const ccomplex_t IT_0223 = (-2)*IT_0032;
    const ccomplex_t IT_0224 = IT_0035 + IT_0223;
    const ccomplex_t IT_0225 = IT_0027*IT_0224;
    const ccomplex_t IT_0226 = (-6)*IT_0225;
    const ccomplex_t IT_0227 = IT_0222 + IT_0226;
    const ccomplex_t IT_0228 = (-18)*IT_0137;
    const ccomplex_t IT_0229 = m_t*IT_0215;
    const ccomplex_t IT_0230 = (-6)*IT_0229;
    const ccomplex_t IT_0231 = pow(s_34, 2);
    const ccomplex_t IT_0232 = IT_0027*IT_0231;
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = IT_0038 + IT_0233;
    const ccomplex_t IT_0235 = m_C_2*IT_0234;
    const ccomplex_t IT_0236 = (-6)*IT_0235;
    const ccomplex_t IT_0237 = conj(IT_0136)*(IT_0068*IT_0181 + IT_0142
      *IT_0186 + IT_0118*IT_0187 + IT_0140*IT_0199 + IT_0163*IT_0210 + IT_0196
      *IT_0221 + IT_0136*IT_0227 + IT_0112*IT_0228 + IT_0169*IT_0230 + IT_0197
      *IT_0236);
    const ccomplex_t IT_0238 = conj(IT_0112)*(IT_0142*IT_0181 + IT_0068
      *IT_0186 + IT_0076*IT_0187 + IT_0140*IT_0200 + IT_0196*IT_0210 + IT_0163
      *IT_0221 + IT_0112*IT_0227 + IT_0136*IT_0228 + IT_0197*IT_0230 + IT_0169
      *IT_0236);
    const ccomplex_t IT_0239 = s_13*IT_0234;
    const ccomplex_t IT_0240 = (-6)*IT_0239;
    const ccomplex_t IT_0241 = IT_0137*IT_0234;
    const ccomplex_t IT_0242 = (-6)*IT_0241;
    const ccomplex_t IT_0243 = (-0.166666666666667)*IT_0068;
    const ccomplex_t IT_0244 = IT_0026*IT_0167 + (-0.166666666666667)*IT_0163
      *IT_0208 + (-0.166666666666667)*IT_0136*IT_0210 + (-0.166666666666667)
      *IT_0196*IT_0212 + (-0.166666666666667)*IT_0169*IT_0217 + (
      -0.166666666666667)*IT_0197*IT_0219 + (-0.166666666666667)*IT_0112*IT_0221
       + IT_0168*IT_0243;
    const ccomplex_t IT_0245 = IT_0026*IT_0173 + (-0.166666666666667)*IT_0163
      *IT_0217 + (-0.166666666666667)*IT_0196*IT_0219 + (-0.166666666666667)
      *IT_0136*IT_0230 + (-0.166666666666667)*IT_0112*IT_0236 + (
      -0.166666666666667)*IT_0169*IT_0240 + (-0.166666666666667)*IT_0197*IT_0242
       + IT_0174*IT_0243;
    const ccomplex_t IT_0246 = conj(IT_0026)*(IT_0026*IT_0037 + IT_0068
      *IT_0070) + (IT_0076*IT_0112 + IT_0118*IT_0136 + IT_0139*IT_0140 + IT_0141
      *IT_0142)*IT_0143 + conj(IT_0068)*(IT_0037*IT_0068 + IT_0026*IT_0070 +
       IT_0163*IT_0168 + IT_0169*IT_0174 + IT_0140*IT_0175 + IT_0142*IT_0176 +
       IT_0136*IT_0181 + IT_0112*IT_0186) + conj(IT_0142)*(IT_0070*IT_0140 +
       IT_0037*IT_0142 + IT_0068*IT_0176 + IT_0112*IT_0181 + IT_0136*IT_0186 +
       IT_0141*IT_0187 + IT_0168*IT_0196 + IT_0174*IT_0197) + (IT_0163*IT_0167 +
       IT_0169*IT_0173)*IT_0198 + conj(IT_0140)*(IT_0037*IT_0140 + IT_0070
      *IT_0142 + IT_0068*IT_0175 + IT_0139*IT_0187 + IT_0136*IT_0199 + IT_0112
      *IT_0200 + IT_0196*IT_0201 + IT_0197*IT_0202) + conj(IT_0196)*(IT_0142
      *IT_0168 + IT_0140*IT_0201 + IT_0196*IT_0208 + IT_0112*IT_0210 + IT_0163
      *IT_0212 + IT_0197*IT_0217 + IT_0169*IT_0219 + IT_0136*IT_0221) + IT_0237 
      + IT_0238 + conj(IT_0197)*(IT_0142*IT_0174 + IT_0140*IT_0202 + IT_0196
      *IT_0217 + IT_0163*IT_0219 + IT_0112*IT_0230 + IT_0136*IT_0236 + IT_0197
      *IT_0240 + IT_0169*IT_0242) + (-6)*conj(IT_0163)*IT_0244 + (-6)*conj
      (IT_0169)*IT_0245;
    return create_ccomplex_return(IT_0246);
}

