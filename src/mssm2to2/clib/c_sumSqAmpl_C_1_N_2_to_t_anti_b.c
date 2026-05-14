#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_2_to_t_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_2_to_t_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = s_13*s_24;
    const ccomplex_t IT_0001 = pow(m_C_1, 2);
    const ccomplex_t IT_0002 = pow(m_N_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0001 + -IT_0002 + -reg_prop, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*(IT_0007 + 0.5*IT_0010);
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0005;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = IT_0003*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = s_34*m_C_1*m_N_2;
    const ccomplex_t IT_0019 = cos(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0025 = IT_0005*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = cos(beta);
    const ccomplex_t IT_0028 = cpow(IT_0027, -1);
    const ccomplex_t IT_0029 = IT_0005*IT_0028;
    const ccomplex_t IT_0030 = pow(m_W, -1);
    const ccomplex_t IT_0031 = m_b*N_d2*e_em*IT_0030*conj(U_sb_10);
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0023 + (-3)*IT_0026 + 3
      *IT_0033);
    const ccomplex_t IT_0035 = 0.166666666666667*IT_0034;
    const ccomplex_t IT_0036 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0037 = IT_0005*IT_0036;
    const ccomplex_t IT_0038 = m_b*conj(U_d1)*V_tb*e_em*IT_0030*U_sb_10;
    const ccomplex_t IT_0039 = IT_0029*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + (-0.5)*IT_0040);
    const ccomplex_t IT_0042 = IT_0035*IT_0041;
    const ccomplex_t IT_0043 = pow(m_t, 2);
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0001 + IT_0043 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = N_B2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0048 = IT_0020*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = N_W2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0051 = IT_0005*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = m_b*N_d2*e_em*IT_0030*conj(U_sb_11);
    const ccomplex_t IT_0054 = IT_0029*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0049 + (-3)*IT_0052 + 3
      *IT_0055);
    const ccomplex_t IT_0057 = 0.166666666666667*IT_0056;
    const ccomplex_t IT_0058 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0059 = IT_0005*IT_0058;
    const ccomplex_t IT_0060 = m_b*conj(U_d1)*V_tb*e_em*IT_0030*U_sb_11;
    const ccomplex_t IT_0061 = IT_0029*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + (-0.5)*IT_0062);
    const ccomplex_t IT_0064 = IT_0057*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0001 + IT_0043 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0046 + -IT_0067;
    const ccomplex_t IT_0069 = m_b*s_13*m_N_2;
    const ccomplex_t IT_0070 = pow(m_W, 2);
    const ccomplex_t IT_0071 = cpow((-2)*s_12 + -IT_0001 + -IT_0002 + IT_0070 
      + -reg_prop, -1);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0005*IT_0030;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = sin(beta);
    const ccomplex_t IT_0075 = IT_0005*IT_0074;
    const ccomplex_t IT_0076 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = 2*IT_0077;
    const ccomplex_t IT_0079 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0080 = IT_0075*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = IT_0020*IT_0074;
    const ccomplex_t IT_0083 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = IT_0081 + IT_0085;
    const ccomplex_t IT_0087 = IT_0078 + IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = IT_0073*IT_0089;
    const ccomplex_t IT_0091 = IT_0071*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = cpow((-2)*s_12 + -IT_0001 + -IT_0002 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0094 = IT_0020*IT_0027;
    const ccomplex_t IT_0095 = IT_0083*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = IT_0005*IT_0027;
    const ccomplex_t IT_0098 = IT_0079*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = IT_0096 + IT_0099;
    const ccomplex_t IT_0101 = IT_0076*IT_0097;
    const ccomplex_t IT_0102 = 2*IT_0101;
    const ccomplex_t IT_0103 = IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0005*IT_0028*IT_0030*IT_0074;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = IT_0105*IT_0107;
    const ccomplex_t IT_0109 = IT_0093*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = IT_0092 + IT_0110;
    const ccomplex_t IT_0112 = m_t*s_14*m_N_2;
    const ccomplex_t IT_0113 = cpow(IT_0074, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0005*IT_0030*IT_0113*U_sb_00;
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = IT_0035*IT_0115;
    const ccomplex_t IT_0117 = IT_0044*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0005*IT_0030*IT_0113*U_sb_01;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = IT_0057*IT_0120;
    const ccomplex_t IT_0122 = IT_0065*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = IT_0118 + IT_0123;
    const ccomplex_t IT_0125 = m_t*s_24*m_C_1;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0005*IT_0030;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0129 = IT_0097*IT_0128;
    const ccomplex_t IT_0130 = (-2)*IT_0129;
    const ccomplex_t IT_0131 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0132 = IT_0094*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0135 = IT_0097*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = IT_0133 + IT_0136;
    const ccomplex_t IT_0138 = IT_0130 + IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = IT_0127*IT_0140;
    const ccomplex_t IT_0142 = IT_0071*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0005*IT_0027*IT_0030*IT_0113;
    const ccomplex_t IT_0145 = 0.5*IT_0144;
    const ccomplex_t IT_0146 = IT_0075*IT_0128;
    const ccomplex_t IT_0147 = (-2)*IT_0146;
    const ccomplex_t IT_0148 = IT_0082*IT_0131;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = IT_0075*IT_0134;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = IT_0149 + IT_0151;
    const ccomplex_t IT_0153 = IT_0147 + IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = 0.5*IT_0154;
    const ccomplex_t IT_0156 = IT_0145*IT_0155;
    const ccomplex_t IT_0157 = IT_0093*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = -IT_0143 + -IT_0158;
    const ccomplex_t IT_0160 = m_b*m_t*s_12;
    const ccomplex_t IT_0161 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0162 = IT_0020*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = m_b*conj(N_d2)*e_em*IT_0030*conj(U_sb_00);
    const ccomplex_t IT_0165 = IT_0029*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0163 + 1.5*IT_0166);
    const ccomplex_t IT_0168 = (-0.333333333333333)*IT_0167;
    const ccomplex_t IT_0169 = IT_0115*IT_0168;
    const ccomplex_t IT_0170 = IT_0044*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = conj(N_B2)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0173 = IT_0020*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = m_b*conj(N_d2)*e_em*IT_0030*conj(U_sb_01);
    const ccomplex_t IT_0176 = IT_0029*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + 1.5*IT_0177);
    const ccomplex_t IT_0179 = (-0.333333333333333)*IT_0178;
    const ccomplex_t IT_0180 = IT_0120*IT_0179;
    const ccomplex_t IT_0181 = IT_0065*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = -IT_0171 + -IT_0182;
    const ccomplex_t IT_0184 = m_b*s_23*m_C_1;
    const ccomplex_t IT_0185 = IT_0041*IT_0168;
    const ccomplex_t IT_0186 = IT_0044*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = IT_0063*IT_0179;
    const ccomplex_t IT_0189 = IT_0065*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = IT_0187 + IT_0190;
    const ccomplex_t IT_0192 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0193 = IT_0005*IT_0192;
    const ccomplex_t IT_0194 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0195 = IT_0005*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + (-0.5)*IT_0196);
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0014*IT_0198;
    const ccomplex_t IT_0200 = IT_0003*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = 12*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0005*IT_0028*IT_0030*conj(U_st_01);
    const ccomplex_t IT_0204 = (-0.5)*IT_0203;
    const ccomplex_t IT_0205 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0206 = IT_0020*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = IT_0005*IT_0113;
    const ccomplex_t IT_0209 = m_t*N_u2*e_em*IT_0030*U_st_01;
    const ccomplex_t IT_0210 = IT_0208*IT_0209;
    const ccomplex_t IT_0211 = 1.4142135623731*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*(IT_0207 + (-0.75)*IT_0211);
    const ccomplex_t IT_0213 = 0.666666666666667*IT_0212;
    const ccomplex_t IT_0214 = IT_0204*IT_0213;
    const ccomplex_t IT_0215 = cpow((-2)*s_23 + IT_0002 + IT_0043 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0216 = IT_0214*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0005*IT_0028*IT_0030*conj(U_st_00);
    const ccomplex_t IT_0219 = (-0.5)*IT_0218;
    const ccomplex_t IT_0220 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0221 = IT_0020*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = m_t*N_u2*e_em*IT_0030*U_st_00;
    const ccomplex_t IT_0224 = IT_0208*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*(IT_0222 + (-0.75)*IT_0225);
    const ccomplex_t IT_0227 = 0.666666666666667*IT_0226;
    const ccomplex_t IT_0228 = IT_0219*IT_0227;
    const ccomplex_t IT_0229 = cpow((-2)*s_23 + IT_0002 + IT_0043 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0230 = IT_0228*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = -IT_0217 + -IT_0231;
    const ccomplex_t IT_0233 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0234 = IT_0020*IT_0233;
    const ccomplex_t IT_0235 = 1.4142135623731*IT_0234;
    const ccomplex_t IT_0236 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0237 = IT_0005*IT_0236;
    const ccomplex_t IT_0238 = 1.4142135623731*IT_0237;
    const ccomplex_t IT_0239 = m_t*conj(N_u2)*e_em*IT_0030*U_st_10;
    const ccomplex_t IT_0240 = IT_0208*IT_0239;
    const ccomplex_t IT_0241 = 1.4142135623731*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*(IT_0235 + 3*IT_0238 + 3
      *IT_0241);
    const ccomplex_t IT_0243 = 0.166666666666667*IT_0242;
    const ccomplex_t IT_0244 = m_t*V_tb*V_u1*e_em*IT_0030*conj(U_st_10);
    const ccomplex_t IT_0245 = IT_0208*IT_0244;
    const ccomplex_t IT_0246 = 1.4142135623731*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = (-0.5)*IT_0247;
    const ccomplex_t IT_0249 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0250 = IT_0005*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = -IT_0252;
    const ccomplex_t IT_0254 = IT_0248 + IT_0253;
    const ccomplex_t IT_0255 = IT_0243*IT_0254;
    const ccomplex_t IT_0256 = IT_0229*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0259 = IT_0020*IT_0258;
    const ccomplex_t IT_0260 = 1.4142135623731*IT_0259;
    const ccomplex_t IT_0261 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0262 = IT_0005*IT_0261;
    const ccomplex_t IT_0263 = 1.4142135623731*IT_0262;
    const ccomplex_t IT_0264 = m_t*conj(N_u2)*e_em*IT_0030*U_st_11;
    const ccomplex_t IT_0265 = IT_0208*IT_0264;
    const ccomplex_t IT_0266 = 1.4142135623731*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*(IT_0260 + 3*IT_0263 + 3
      *IT_0266);
    const ccomplex_t IT_0268 = 0.166666666666667*IT_0267;
    const ccomplex_t IT_0269 = m_t*V_tb*V_u1*e_em*IT_0030*conj(U_st_11);
    const ccomplex_t IT_0270 = IT_0208*IT_0269;
    const ccomplex_t IT_0271 = 1.4142135623731*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = (-0.5)*IT_0272;
    const ccomplex_t IT_0274 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0275 = IT_0005*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*IT_0275;
    const ccomplex_t IT_0277 = -IT_0276;
    const ccomplex_t IT_0278 = -IT_0277;
    const ccomplex_t IT_0279 = IT_0273 + IT_0278;
    const ccomplex_t IT_0280 = IT_0268*IT_0279;
    const ccomplex_t IT_0281 = IT_0215*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = IT_0257 + IT_0282;
    const ccomplex_t IT_0284 = IT_0089*IT_0127;
    const ccomplex_t IT_0285 = IT_0071*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = IT_0105*IT_0145;
    const ccomplex_t IT_0288 = IT_0093*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = IT_0286 + IT_0289;
    const ccomplex_t IT_0291 = IT_0227*IT_0254;
    const ccomplex_t IT_0292 = IT_0229*IT_0291;
    const ccomplex_t IT_0293 = (0 + _Complex_I*1)*IT_0292;
    const ccomplex_t IT_0294 = IT_0213*IT_0279;
    const ccomplex_t IT_0295 = IT_0215*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = -IT_0293 + -IT_0296;
    const ccomplex_t IT_0298 = conj(IT_0290) + conj(IT_0297);
    const ccomplex_t IT_0299 = 12*IT_0017;
    const ccomplex_t IT_0300 = IT_0204*IT_0268;
    const ccomplex_t IT_0301 = IT_0215*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*IT_0301;
    const ccomplex_t IT_0303 = IT_0219*IT_0243;
    const ccomplex_t IT_0304 = IT_0229*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*IT_0304;
    const ccomplex_t IT_0306 = IT_0302 + IT_0305;
    const ccomplex_t IT_0307 = IT_0069*IT_0299;
    const ccomplex_t IT_0308 = m_b*m_t*m_C_1*m_N_2;
    const ccomplex_t IT_0309 = s_14*s_23;
    const ccomplex_t IT_0310 = 24*IT_0201;
    const ccomplex_t IT_0311 = 24*IT_0017;
    const ccomplex_t IT_0312 = 24*conj(IT_0017);
    const ccomplex_t IT_0313 = IT_0073*IT_0140;
    const ccomplex_t IT_0314 = IT_0071*IT_0313;
    const ccomplex_t IT_0315 = (0 + _Complex_I*1)*IT_0314;
    const ccomplex_t IT_0316 = IT_0107*IT_0155;
    const ccomplex_t IT_0317 = IT_0093*IT_0316;
    const ccomplex_t IT_0318 = (0 + _Complex_I*1)*IT_0317;
    const ccomplex_t IT_0319 = -IT_0315 + -IT_0318;
    const ccomplex_t IT_0320 = (-12)*IT_0201;
    const ccomplex_t IT_0321 = (-12)*IT_0017;
    const ccomplex_t IT_0322 = s_12*s_34;
    const ccomplex_t IT_0323 = -IT_0309;
    const ccomplex_t IT_0324 = IT_0322 + IT_0323;
    const ccomplex_t IT_0325 = IT_0000 + IT_0324;
    const ccomplex_t IT_0326 = 6*IT_0325;
    const ccomplex_t IT_0327 = 6*IT_0184;
    const ccomplex_t IT_0328 = 6*IT_0308;
    const ccomplex_t IT_0329 = 6*IT_0112;
    const ccomplex_t IT_0330 = 12*IT_0000;
    const ccomplex_t IT_0331 = 12*IT_0125;
    const ccomplex_t IT_0332 = 12*IT_0308;
    const ccomplex_t IT_0333 = 12*IT_0069;
    const ccomplex_t IT_0334 = 12*conj(IT_0201);
    const ccomplex_t IT_0335 = (-6)*IT_0160;
    const ccomplex_t IT_0336 = -IT_0000;
    const ccomplex_t IT_0337 = IT_0324 + IT_0336;
    const ccomplex_t IT_0338 = (-6)*IT_0337;
    const ccomplex_t IT_0339 = (-6)*IT_0018;
    const ccomplex_t IT_0340 = IT_0069*IT_0312 + conj(IT_0319)*IT_0326 + conj
      (IT_0283)*IT_0327 + IT_0298*IT_0328 + conj(IT_0232)*IT_0329 + conj(IT_0191
      )*IT_0330 + conj(IT_0183)*IT_0331 + conj(IT_0124)*IT_0332 + conj(IT_0068)
      *IT_0333 + IT_0184*IT_0334 + conj(IT_0159)*IT_0335 + conj(IT_0306)*IT_0338
       + conj(IT_0111)*IT_0339;
    const ccomplex_t IT_0341 = conj(IT_0306) + conj(IT_0319);
    const ccomplex_t IT_0342 = IT_0125*IT_0312 + conj(IT_0290)*IT_0326 + conj
      (IT_0232)*IT_0327 + conj(IT_0283)*IT_0329 + conj(IT_0124)*IT_0330 + conj
      (IT_0068)*IT_0331 + conj(IT_0191)*IT_0332 + conj(IT_0183)*IT_0333 +
       IT_0112*IT_0334 + conj(IT_0111)*IT_0335 + conj(IT_0297)*IT_0338 + conj
      (IT_0159)*IT_0339 + IT_0328*IT_0341;
    const ccomplex_t IT_0343 = 6*IT_0160;
    const ccomplex_t IT_0344 = 6*IT_0125;
    const ccomplex_t IT_0345 = 6*IT_0018;
    const ccomplex_t IT_0346 = 6*IT_0069;
    const ccomplex_t IT_0347 = 12*IT_0309;
    const ccomplex_t IT_0348 = 12*IT_0184;
    const ccomplex_t IT_0349 = 12*IT_0112;
    const ccomplex_t IT_0350 = 12*conj(IT_0017);
    const ccomplex_t IT_0351 = 24*conj(IT_0201);
    const ccomplex_t IT_0352 = (-6)*IT_0184;
    const ccomplex_t IT_0353 = (-6)*IT_0112;
    const ccomplex_t IT_0354 = conj(IT_0191)*IT_0327 + conj(IT_0124)*IT_0329 +
       conj(IT_0232)*IT_0332 + conj(IT_0183)*IT_0343 + conj(IT_0159)*IT_0344 +
       conj(IT_0068)*IT_0345 + conj(IT_0111)*IT_0346 + conj(IT_0283)*IT_0347 +
       conj(IT_0306)*IT_0348 + conj(IT_0297)*IT_0349 + IT_0018*IT_0350 + IT_0309
      *IT_0351 + conj(IT_0319)*IT_0352 + conj(IT_0290)*IT_0353;
    const ccomplex_t IT_0355 = IT_0308*IT_0312 + conj(IT_0183)*IT_0330 + conj
      (IT_0191)*IT_0331 + conj(IT_0068)*IT_0332 + conj(IT_0124)*IT_0333 +
       IT_0160*IT_0334 + conj(IT_0283)*IT_0343 + IT_0341*IT_0344 + conj(IT_0232)
      *IT_0345 + IT_0298*IT_0346 + conj(IT_0159)*IT_0352 + conj(IT_0111)*IT_0353;
    const ccomplex_t IT_0356 = IT_0000*IT_0312 + conj(IT_0068)*IT_0330 + conj
      (IT_0124)*IT_0331 + conj(IT_0183)*IT_0332 + conj(IT_0191)*IT_0333 +
       IT_0018*IT_0334 + conj(IT_0232)*IT_0343 + IT_0298*IT_0344 + conj(IT_0283)
      *IT_0345 + IT_0341*IT_0346 + conj(IT_0111)*IT_0352 + conj(IT_0159)*IT_0353;
    const ccomplex_t IT_0357 = conj(IT_0124)*IT_0327 + conj(IT_0191)*IT_0329 +
       conj(IT_0283)*IT_0332 + conj(IT_0068)*IT_0343 + conj(IT_0111)*IT_0344 +
       conj(IT_0183)*IT_0345 + conj(IT_0159)*IT_0346 + conj(IT_0232)*IT_0347 +
       conj(IT_0297)*IT_0348 + conj(IT_0306)*IT_0349 + IT_0160*IT_0350 + IT_0308
      *IT_0351 + conj(IT_0290)*IT_0352 + conj(IT_0319)*IT_0353;
    const ccomplex_t IT_0358 = (-12)*IT_0160;
    const ccomplex_t IT_0359 = (-12)*IT_0018;
    const ccomplex_t IT_0360 = (-12)*conj(IT_0017);
    const ccomplex_t IT_0361 = 12*IT_0322;
    const ccomplex_t IT_0362 = conj(IT_0159)*IT_0332 + IT_0069*IT_0334 + conj
      (IT_0124)*IT_0335 + conj(IT_0191)*IT_0339 + conj(IT_0297)*IT_0343 + conj
      (IT_0232)*IT_0344 + conj(IT_0306)*IT_0345 + conj(IT_0283)*IT_0346 + conj
      (IT_0068)*IT_0352 + conj(IT_0183)*IT_0353 + conj(IT_0290)*IT_0358 + conj
      (IT_0319)*IT_0359 + IT_0184*IT_0360 + conj(IT_0111)*IT_0361;
    const ccomplex_t IT_0363 = conj(IT_0111)*IT_0332 + IT_0125*IT_0334 + conj
      (IT_0191)*IT_0335 + conj(IT_0124)*IT_0339 + conj(IT_0306)*IT_0343 + conj
      (IT_0283)*IT_0344 + conj(IT_0297)*IT_0345 + conj(IT_0232)*IT_0346 + conj
      (IT_0183)*IT_0352 + conj(IT_0068)*IT_0353 + conj(IT_0319)*IT_0358 + conj
      (IT_0290)*IT_0359 + IT_0112*IT_0360 + conj(IT_0159)*IT_0361;
    const ccomplex_t IT_0364 = conj(IT_0124)*IT_0328;
    const ccomplex_t IT_0365 = conj(IT_0068)*IT_0346;
    const ccomplex_t IT_0366 = IT_0069*IT_0350;
    const ccomplex_t IT_0367 = IT_0309 + IT_0322 + IT_0336;
    const ccomplex_t IT_0368 = (-6)*IT_0367;
    const ccomplex_t IT_0369 = (-6)*IT_0308;
    const ccomplex_t IT_0370 = conj(IT_0183)*IT_0344;
    const ccomplex_t IT_0371 = conj(IT_0297)*IT_0332 + conj(IT_0191)*IT_0338 +
       conj(IT_0159)*IT_0343 + conj(IT_0111)*IT_0345 + conj(IT_0306)*IT_0347 +
       conj(IT_0283)*IT_0348 + conj(IT_0232)*IT_0349 + IT_0184*IT_0351 + IT_0364
       + IT_0365 + IT_0366 + conj(IT_0319)*IT_0368 + conj(IT_0290)*IT_0369 +
       IT_0370;
    const ccomplex_t IT_0372 = (-12)*conj(IT_0201);
    const ccomplex_t IT_0373 = conj(IT_0191)*IT_0326 + conj(IT_0290)*IT_0332 +
       conj(IT_0283)*IT_0352 + conj(IT_0232)*IT_0353 + conj(IT_0159)*IT_0358 +
       conj(IT_0111)*IT_0359 + conj(IT_0319)*IT_0361 + IT_0364 + IT_0365 +
       IT_0366 + conj(IT_0306)*IT_0368 + conj(IT_0297)*IT_0369 + IT_0370 +
       IT_0184*IT_0372;
    const ccomplex_t IT_0374 = IT_0125*IT_0350;
    const ccomplex_t IT_0375 = conj(IT_0068)*IT_0344;
    const ccomplex_t IT_0376 = conj(IT_0183)*IT_0346;
    const ccomplex_t IT_0377 = conj(IT_0191)*IT_0328;
    const ccomplex_t IT_0378 = conj(IT_0306)*IT_0332 + conj(IT_0124)*IT_0338 +
       conj(IT_0111)*IT_0343 + conj(IT_0159)*IT_0345 + conj(IT_0297)*IT_0347 +
       conj(IT_0232)*IT_0348 + conj(IT_0283)*IT_0349 + IT_0112*IT_0351 + conj
      (IT_0290)*IT_0368 + conj(IT_0319)*IT_0369 + IT_0374 + IT_0375 + IT_0376 +
       IT_0377;
    const ccomplex_t IT_0379 = conj(IT_0124)*IT_0326 + conj(IT_0319)*IT_0332 +
       conj(IT_0232)*IT_0352 + conj(IT_0283)*IT_0353 + conj(IT_0111)*IT_0358 +
       conj(IT_0159)*IT_0359 + conj(IT_0290)*IT_0361 + conj(IT_0297)*IT_0368 +
       conj(IT_0306)*IT_0369 + IT_0112*IT_0372 + IT_0374 + IT_0375 + IT_0376 +
       IT_0377;
    const ccomplex_t IT_0380 = 48*IT_0000*IT_0017*conj(IT_0017) + (IT_0018
      *conj(IT_0068) + IT_0069*conj(IT_0111) + IT_0112*conj(IT_0124) + IT_0125
      *conj(IT_0159) + IT_0160*conj(IT_0183) + IT_0184*conj(IT_0191))*IT_0202 + 
      (IT_0160*conj(IT_0232) + IT_0018*conj(IT_0283) + IT_0125*IT_0298)*IT_0299 
      + conj(IT_0306)*IT_0307 + (IT_0112*conj(IT_0297) + IT_0184*conj(IT_0306) +
       conj(IT_0232)*IT_0308 + conj(IT_0283)*IT_0309)*IT_0310 + (IT_0000*conj
      (IT_0068) + conj(IT_0124)*IT_0125 + IT_0069*conj(IT_0191) + IT_0018*conj
      (IT_0201) + conj(IT_0183)*IT_0308)*IT_0311 + IT_0201*(48*conj(IT_0201)
      *IT_0309 + IT_0018*IT_0312) + conj(IT_0319)*(IT_0307 + IT_0184*IT_0320) +
       conj(IT_0111)*IT_0184*IT_0321 + IT_0112*(conj(IT_0290)*IT_0320 + conj
      (IT_0159)*IT_0321) + IT_0191*IT_0340 + IT_0124*IT_0342 + IT_0283*IT_0354 +
       IT_0183*IT_0355 + IT_0068*IT_0356 + IT_0232*IT_0357 + IT_0111*IT_0362 +
       IT_0159*IT_0363 + IT_0306*IT_0371 + IT_0319*IT_0373 + IT_0297*IT_0378 +
       IT_0290*IT_0379;
    return create_ccomplex_return(IT_0380);
}

