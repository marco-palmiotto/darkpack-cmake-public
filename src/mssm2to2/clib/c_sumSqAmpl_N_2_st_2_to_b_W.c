#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_st_2_to_b_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_st_2_to_b_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
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
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
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
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(IT_0003 + 0.5*IT_0006);
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0001*IT_0009*IT_0011*U_st_01;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    const ccomplex_t IT_0015 = pow(m_b, 2);
    const ccomplex_t IT_0016 = pow(m_st_2, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0021 = IT_0001*IT_0020;
    const ccomplex_t IT_0022 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0023 = IT_0001*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0021 + 0.5*IT_0024);
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0001*IT_0009*IT_0011*U_st_01;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_00)*U_st_01;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = cos(theta_W);
    const ccomplex_t IT_0036 = cpow(IT_0035, -1);
    const ccomplex_t IT_0037 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = IT_0001*IT_0011;
    const ccomplex_t IT_0041 = m_b*N_d2*e_em*IT_0009*U_sb_00;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0039 + 1.5*IT_0043);
    const ccomplex_t IT_0045 = 0.666666666666667*IT_0044;
    const ccomplex_t IT_0046 = 1.5*IT_0045;
    const ccomplex_t IT_0047 = 0.666666666666667*IT_0046;
    const ccomplex_t IT_0048 = 1.5*IT_0047;
    const ccomplex_t IT_0049 = 0.666666666666667*IT_0048;
    const ccomplex_t IT_0050 = 1.5*IT_0049;
    const ccomplex_t IT_0051 = 0.666666666666667*IT_0050;
    const ccomplex_t IT_0052 = 1.5*IT_0051;
    const ccomplex_t IT_0053 = 0.666666666666667*IT_0052;
    const ccomplex_t IT_0054 = 1.5*IT_0053;
    const ccomplex_t IT_0055 = 0.666666666666667*IT_0054;
    const ccomplex_t IT_0056 = 1.5*IT_0055;
    const ccomplex_t IT_0057 = (-0.333333333333333)*IT_0056;
    const ccomplex_t IT_0058 = pow(m_N_2, 2);
    const ccomplex_t IT_0059 = cpow((-2)*s_13 + IT_0015 + IT_0058 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = IT_0034*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_01)*U_st_01;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = m_b*N_d2*e_em*IT_0009*U_sb_01;
    const ccomplex_t IT_0066 = IT_0040*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0069 = IT_0036*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0067 + 0.666666666666667
      *IT_0070);
    const ccomplex_t IT_0072 = 1.5*IT_0071;
    const ccomplex_t IT_0073 = (-0.333333333333333)*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_13 + IT_0015 + IT_0058 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0064*IT_0076;
    const ccomplex_t IT_0078 = (-2)*IT_0019 + (-2)*IT_0032 + IT_0062 + IT_0077;
    const ccomplex_t IT_0079 = pow(m_W, -2);
    const ccomplex_t IT_0080 = pow(s_14, 2);
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = -IT_0058;
    const ccomplex_t IT_0083 = IT_0081 + IT_0082;
    const ccomplex_t IT_0084 = s_13*IT_0083;
    const ccomplex_t IT_0085 = m_b*conj(N_d2)*e_em*IT_0009*U_sb_11;
    const ccomplex_t IT_0086 = IT_0040*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0089 = IT_0001*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0092 = IT_0036*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0087 + -IT_0090 +
       0.333333333333333*IT_0093);
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = (-3)*IT_0095;
    const ccomplex_t IT_0097 = 0.166666666666667*IT_0096;
    const ccomplex_t IT_0098 = IT_0074*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = IT_0064*IT_0099;
    const ccomplex_t IT_0101 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0102 = IT_0036*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = m_b*conj(N_d2)*e_em*IT_0009*U_sb_10;
    const ccomplex_t IT_0105 = IT_0040*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0108 = IT_0001*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + 3*IT_0106 + (-3)
      *IT_0109);
    const ccomplex_t IT_0111 = 0.333333333333333*IT_0110;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = 3*IT_0115;
    const ccomplex_t IT_0117 = 0.333333333333333*IT_0116;
    const ccomplex_t IT_0118 = 3*IT_0117;
    const ccomplex_t IT_0119 = (-0.333333333333333)*IT_0118;
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = 3*IT_0120;
    const ccomplex_t IT_0122 = 0.166666666666667*IT_0121;
    const ccomplex_t IT_0123 = IT_0059*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = IT_0034*IT_0124;
    const ccomplex_t IT_0126 = -IT_0100 + -IT_0125;
    const ccomplex_t IT_0127 = m_b*m_N_2;
    const ccomplex_t IT_0128 = IT_0083*IT_0127;
    const ccomplex_t IT_0129 = IT_0019 + IT_0032;
    const ccomplex_t IT_0130 = s_13*s_14;
    const ccomplex_t IT_0131 = s_34*IT_0058;
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = IT_0130 + IT_0132;
    const ccomplex_t IT_0134 = cpow((-2)*s_12 + -IT_0016 + -IT_0058 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0135 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0136 = IT_0036*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0139 = IT_0001*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = sin(beta);
    const ccomplex_t IT_0142 = cpow(IT_0141, -1);
    const ccomplex_t IT_0143 = IT_0001*IT_0142;
    const ccomplex_t IT_0144 = m_t*conj(N_u2)*e_em*IT_0009*U_st_11;
    const ccomplex_t IT_0145 = IT_0143*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0137 + 3*IT_0140 + 3
      *IT_0146);
    const ccomplex_t IT_0148 = 0.166666666666667*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = IT_0148*IT_0150;
    const ccomplex_t IT_0152 = IT_0134*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = pow(m_W, 2);
    const ccomplex_t IT_0155 = IT_0127*IT_0154;
    const ccomplex_t IT_0156 = 3*IT_0155;
    const ccomplex_t IT_0157 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0158 = IT_0001*IT_0157;
    const ccomplex_t IT_0159 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0160 = IT_0001*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*(IT_0158 + (-0.5)*IT_0161);
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = V_tb*e_em*conj(V_Wp2)*U_st_01;
    const ccomplex_t IT_0165 = IT_0001*IT_0164;
    const ccomplex_t IT_0166 = m_t*V_tb*conj(V_u2)*e_em*IT_0009*U_st_11;
    const ccomplex_t IT_0167 = IT_0143*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(IT_0165 + (-0.5)*IT_0168);
    const ccomplex_t IT_0170 = -IT_0169;
    const ccomplex_t IT_0171 = IT_0163*IT_0170;
    const ccomplex_t IT_0172 = IT_0030*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0175 = IT_0001*IT_0174;
    const ccomplex_t IT_0176 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0177 = IT_0001*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*(IT_0175 + (-0.5)*IT_0178);
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = V_tb*e_em*conj(V_Wp1)*U_st_01;
    const ccomplex_t IT_0182 = IT_0001*IT_0181;
    const ccomplex_t IT_0183 = m_t*V_tb*conj(V_u1)*e_em*IT_0009*U_st_11;
    const ccomplex_t IT_0184 = IT_0143*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0182 + (-0.5)*IT_0185);
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = IT_0180*IT_0187;
    const ccomplex_t IT_0189 = IT_0017*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = -IT_0173 + -IT_0190;
    const ccomplex_t IT_0192 = (-3)*IT_0155;
    const ccomplex_t IT_0193 = -IT_0062 + -IT_0077;
    const ccomplex_t IT_0194 = s_14*IT_0015;
    const ccomplex_t IT_0195 = s_13*s_34;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = pow(m_W, 4);
    const ccomplex_t IT_0199 = s_13*IT_0198;
    const ccomplex_t IT_0200 = IT_0079*IT_0199;
    const ccomplex_t IT_0201 = s_14*s_34;
    const ccomplex_t IT_0202 = 24*IT_0201;
    const ccomplex_t IT_0203 = s_13*IT_0154;
    const ccomplex_t IT_0204 = (-12)*IT_0203;
    const ccomplex_t IT_0205 = 6*IT_0200 + IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = s_14*IT_0154;
    const ccomplex_t IT_0207 = IT_0079*IT_0206;
    const ccomplex_t IT_0208 = 0.5*IT_0207;
    const ccomplex_t IT_0209 = s_14 + IT_0208;
    const ccomplex_t IT_0210 = m_b*IT_0209;
    const ccomplex_t IT_0211 = 2*IT_0210;
    const ccomplex_t IT_0212 = IT_0008*IT_0187;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0214 = IT_0017*IT_0212*IT_0213;
    const ccomplex_t IT_0215 = IT_0026*IT_0170;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0217 = IT_0030*IT_0215*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0219 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0220 = IT_0036*IT_0219;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = m_t*N_u2*e_em*IT_0009*U_st_01;
    const ccomplex_t IT_0223 = IT_0143*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*(IT_0221 + (-0.75)*IT_0224);
    const ccomplex_t IT_0226 = 0.666666666666667*IT_0225;
    const ccomplex_t IT_0227 = IT_0150*IT_0226;
    const ccomplex_t IT_0228 = IT_0134*IT_0218*IT_0227;
    const ccomplex_t IT_0229 = -IT_0214 + -IT_0217 + -IT_0228;
    const ccomplex_t IT_0230 = IT_0173 + IT_0190;
    const ccomplex_t IT_0231 = m_N_2*IT_0230;
    const ccomplex_t IT_0232 = IT_0229 + IT_0231;
    const ccomplex_t IT_0233 = s_34*IT_0154;
    const ccomplex_t IT_0234 = IT_0079*IT_0233;
    const ccomplex_t IT_0235 = (-0.25)*IT_0234;
    const ccomplex_t IT_0236 = s_34 + IT_0235;
    const ccomplex_t IT_0237 = m_N_2*IT_0236;
    const ccomplex_t IT_0238 = (-4)*IT_0237;
    const ccomplex_t IT_0239 = m_N_2*IT_0032;
    const ccomplex_t IT_0240 = m_N_2*IT_0019;
    const ccomplex_t IT_0241 = m_b*IT_0153;
    const ccomplex_t IT_0242 = IT_0239 + IT_0240 + -IT_0241;
    const ccomplex_t IT_0243 = IT_0028*IT_0163;
    const ccomplex_t IT_0244 = IT_0030*IT_0216*IT_0243;
    const ccomplex_t IT_0245 = IT_0013*IT_0180;
    const ccomplex_t IT_0246 = IT_0017*IT_0213*IT_0245;
    const ccomplex_t IT_0247 = -IT_0244 + -IT_0246;
    const ccomplex_t IT_0248 = IT_0242 + IT_0247;
    const ccomplex_t IT_0249 = (-18)*conj(IT_0229);
    const ccomplex_t IT_0250 = s_14*s_34*IT_0079;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = s_13 + IT_0251;
    const ccomplex_t IT_0253 = s_13*IT_0252;
    const ccomplex_t IT_0254 = 2*IT_0173 + 2*IT_0190;
    const ccomplex_t IT_0255 = IT_0100 + IT_0125 + 2*IT_0153;
    const ccomplex_t IT_0256 = IT_0127*IT_0252;
    const ccomplex_t IT_0257 = m_N_2*IT_0252;
    const ccomplex_t IT_0258 = (-6)*conj(IT_0078);
    const ccomplex_t IT_0259 = 0.5*IT_0234;
    const ccomplex_t IT_0260 = s_34 + IT_0259;
    const ccomplex_t IT_0261 = m_N_2*IT_0260;
    const ccomplex_t IT_0262 = 2*IT_0261;
    const ccomplex_t IT_0263 = (-2)*IT_0201;
    const ccomplex_t IT_0264 = IT_0203 + IT_0263;
    const ccomplex_t IT_0265 = IT_0079*IT_0264;
    const ccomplex_t IT_0266 = (-6)*IT_0265;
    const ccomplex_t IT_0267 = 12*s_13;
    const ccomplex_t IT_0268 = IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = m_b*IT_0252;
    const ccomplex_t IT_0270 = (-6)*IT_0193;
    const ccomplex_t IT_0271 = pow(s_34, 2);
    const ccomplex_t IT_0272 = IT_0079*IT_0271;
    const ccomplex_t IT_0273 = -IT_0272;
    const ccomplex_t IT_0274 = IT_0015 + IT_0273;
    const ccomplex_t IT_0275 = IT_0127*IT_0274;
    const ccomplex_t IT_0276 = m_N_2*IT_0274;
    const ccomplex_t IT_0277 = (-6)*conj(IT_0193);
    const ccomplex_t IT_0278 = (-6)*IT_0200;
    const ccomplex_t IT_0279 = IT_0202 + IT_0278;
    const ccomplex_t IT_0280 = (-6)*IT_0126;
    const ccomplex_t IT_0281 = 0.166666666666667*IT_0280;
    const ccomplex_t IT_0282 = IT_0254 + IT_0281;
    const ccomplex_t IT_0283 = (-0.25)*IT_0207;
    const ccomplex_t IT_0284 = s_14 + IT_0283;
    const ccomplex_t IT_0285 = m_b*IT_0284;
    const ccomplex_t IT_0286 = (-4)*IT_0285;
    const ccomplex_t IT_0287 = (-6)*IT_0255;
    const ccomplex_t IT_0288 = IT_0153 + IT_0191;
    const ccomplex_t IT_0289 = (-6)*conj(IT_0255);
    const ccomplex_t IT_0290 = (-6)*IT_0229;
    const ccomplex_t IT_0291 = (-6)*conj(IT_0191)*(IT_0129*IT_0192 + (
      -0.166666666666667)*IT_0191*IT_0205 + IT_0211*IT_0248 + (
      -0.166666666666667)*IT_0153*IT_0279 + -IT_0133*IT_0282 + (
      -0.166666666666667)*IT_0197*IT_0287 + IT_0238*(IT_0231 + (
      -0.166666666666667)*IT_0290));
    const ccomplex_t IT_0292 = (-6)*conj(IT_0229);
    const ccomplex_t IT_0293 = (-6)*conj(IT_0254);
    const ccomplex_t IT_0294 = m_b*IT_0083;
    const ccomplex_t IT_0295 = 6*IT_0294;
    const ccomplex_t IT_0296 = (-6)*IT_0078;
    const ccomplex_t IT_0297 = 6*IT_0129;
    const ccomplex_t IT_0298 = (-6)*IT_0254;
    const ccomplex_t IT_0299 = 6*IT_0153;
    const ccomplex_t IT_0300 = (-6)*IT_0191;
    const ccomplex_t IT_0301 = IT_0127*((-18)*IT_0229 + (-18)*IT_0231) +
       IT_0248*IT_0268 + IT_0270*IT_0276 + IT_0269*IT_0287 + IT_0257*IT_0296 +
       IT_0238*IT_0297 + IT_0295*(IT_0126 + 0.166666666666667*IT_0298) + IT_0286
      *IT_0299 + IT_0211*IT_0300;
    const ccomplex_t IT_0302 = IT_0127*((-18)*IT_0242 + (-18)*IT_0247) +
       IT_0232*IT_0268 + IT_0269*IT_0270 + 6*IT_0257*IT_0282 + IT_0276*IT_0287 +
       6*IT_0078*IT_0294 + IT_0211*IT_0297 + IT_0262*IT_0299 + IT_0238*IT_0300;
    const ccomplex_t IT_0303 = 0.166666666666667*IT_0256;
    const ccomplex_t IT_0304 = s_13*IT_0274;
    const ccomplex_t IT_0305 = 0.166666666666667*IT_0304;
    const ccomplex_t IT_0306 = (-6)*conj(IT_0126);
    const ccomplex_t IT_0307 = IT_0133*IT_0288;
    const ccomplex_t IT_0308 = 6*conj(IT_0078)*(IT_0078*IT_0084 + IT_0126
      *IT_0128 + IT_0129*IT_0133) + 6*conj(IT_0129)*(IT_0078*IT_0133 + IT_0153
      *IT_0156 + -IT_0191*IT_0192 + IT_0193*IT_0197 + 0.166666666666667*IT_0129
      *IT_0205 + IT_0211*IT_0232 + IT_0238*IT_0248) + IT_0127*IT_0248*IT_0249 + 
      (IT_0193*IT_0253 + IT_0128*IT_0254 + IT_0255*IT_0256 + IT_0248*IT_0257)
      *IT_0258 + 6*conj(IT_0229)*(IT_0129*IT_0211 + IT_0254*IT_0257 + IT_0153
      *IT_0262 + 0.166666666666667*IT_0232*IT_0268 + 0.166666666666667*IT_0269
      *IT_0270) + (IT_0126*IT_0256 + IT_0232*IT_0269 + IT_0255*IT_0275 + IT_0248
      *IT_0276)*IT_0277 + 6*conj(IT_0153)*(IT_0129*IT_0156 + 0.166666666666667
      *IT_0153*IT_0205 + IT_0232*IT_0262 + 0.166666666666667*IT_0191*IT_0279 +
       IT_0133*IT_0282 + IT_0248*IT_0286 + 0.166666666666667*IT_0197*IT_0287) + 
      (IT_0248*IT_0269 + IT_0231*IT_0276 + IT_0197*IT_0288)*IT_0289 + IT_0291 + 
      (IT_0191*IT_0238 + IT_0126*IT_0257 + IT_0255*IT_0276)*IT_0292 + (IT_0078
      *conj(IT_0229) + conj(IT_0078)*IT_0232 + IT_0248*(conj(IT_0126) +
       0.166666666666667*IT_0293))*IT_0295 + (conj(IT_0242) + conj(IT_0247))
      *IT_0301 + conj(IT_0231)*IT_0302 + 6*conj(IT_0126)*(IT_0084*IT_0126 +
       IT_0078*IT_0128 + 0.166666666666667*IT_0257*IT_0290 + IT_0270*IT_0303) +
       6*conj(IT_0193)*(IT_0129*IT_0197 + IT_0254*IT_0256 + 0.166666666666667
      *IT_0253*IT_0296 + IT_0270*IT_0305) + 6*conj(IT_0255)*(0.166666666666667
      *IT_0270*IT_0275 + IT_0253*IT_0282 + 0.166666666666667*IT_0276*IT_0290 +
       IT_0296*IT_0303 + IT_0287*IT_0305) + IT_0306*(IT_0084*IT_0254 + IT_0253
      *IT_0255 + IT_0231*IT_0257 + IT_0307) + 6*conj(IT_0254)*(IT_0253*IT_0255 +
       IT_0193*IT_0256 + IT_0232*IT_0257 + IT_0084*IT_0282 + 0.166666666666667
      *IT_0128*IT_0296 + IT_0307);
    return create_ccomplex_return(IT_0308);
}

