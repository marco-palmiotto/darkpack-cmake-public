#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_sb_1_to_t_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_sb_1_to_t_anti_W(
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
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
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
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_t, 2);
    const ccomplex_t IT_0001 = pow(m_sb_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + (-0.5)*IT_0009);
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0004*IT_0012*IT_0014*U_sb_00;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0011*IT_0016;
    const ccomplex_t IT_0018 = IT_0002*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0021 = N_W4*e_em*conj(V_Wp1);
    const ccomplex_t IT_0022 = IT_0004*IT_0021;
    const ccomplex_t IT_0023 = N_u4*conj(V_u1)*e_em;
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0022 + (-0.5)*IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0004*IT_0012*IT_0014*U_sb_00;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*U_sb_00*conj(U_st_01);
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = cos(theta_W);
    const ccomplex_t IT_0036 = cpow(IT_0035, -1);
    const ccomplex_t IT_0037 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = IT_0004*IT_0014;
    const ccomplex_t IT_0041 = m_t*N_u4*e_em*IT_0012*U_st_01;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0039 + (-0.75)*IT_0043);
    const ccomplex_t IT_0045 = 0.666666666666667*IT_0044;
    const ccomplex_t IT_0046 = pow(m_N_4, 2);
    const ccomplex_t IT_0047 = cpow((-2)*s_13 + IT_0000 + IT_0046 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0034*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*U_sb_00*conj(U_st_00);
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0054 = IT_0036*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = m_t*N_u4*e_em*IT_0012*U_st_00;
    const ccomplex_t IT_0057 = IT_0040*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + (-0.75)*IT_0058);
    const ccomplex_t IT_0060 = 0.666666666666667*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_13 + IT_0000 + IT_0046 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0063;
    const ccomplex_t IT_0065 = 2*IT_0019 + 2*IT_0032 + IT_0050 + IT_0064;
    const ccomplex_t IT_0066 = pow(m_W, -2);
    const ccomplex_t IT_0067 = pow(s_14, 2);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = -IT_0046;
    const ccomplex_t IT_0070 = IT_0068 + IT_0069;
    const ccomplex_t IT_0071 = s_13*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0073 = cos(beta);
    const ccomplex_t IT_0074 = cpow(IT_0073, -1);
    const ccomplex_t IT_0075 = IT_0004*IT_0074;
    const ccomplex_t IT_0076 = m_b*conj(U_d2)*V_tb*e_em*IT_0012*U_sb_10;
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0080 = IT_0004*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0078 + (-2)*IT_0080);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0011*IT_0082;
    const ccomplex_t IT_0084 = IT_0002*IT_0072*IT_0083;
    const ccomplex_t IT_0085 = pow(m_b, 2);
    const ccomplex_t IT_0086 = cpow((-2)*s_12 + -IT_0001 + -IT_0046 + IT_0085 
      + -reg_prop, -1);
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0091 = IT_0036*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = m_b*N_d4*e_em*IT_0012*U_sb_00;
    const ccomplex_t IT_0094 = IT_0075*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0092 + 1.5*IT_0095);
    const ccomplex_t IT_0097 = (-0.333333333333333)*IT_0096;
    const ccomplex_t IT_0098 = IT_0089*IT_0097;
    const ccomplex_t IT_0099 = IT_0086*IT_0087*IT_0098;
    const ccomplex_t IT_0100 = conj(N_W4)*e_em*U_Wm1;
    const ccomplex_t IT_0101 = IT_0004*IT_0100;
    const ccomplex_t IT_0102 = conj(N_d4)*U_d1*e_em;
    const ccomplex_t IT_0103 = IT_0004*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + 0.5*IT_0104);
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = m_b*conj(U_d1)*V_tb*e_em*IT_0012*U_sb_10;
    const ccomplex_t IT_0108 = IT_0075*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0111 = IT_0004*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0109 + (-2)*IT_0111);
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0106*IT_0113;
    const ccomplex_t IT_0115 = IT_0020*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = m_N_4*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0119 = IT_0027*IT_0113;
    const ccomplex_t IT_0120 = IT_0020*IT_0118*IT_0119;
    const ccomplex_t IT_0121 = -IT_0084 + -IT_0099 + IT_0117 + -IT_0120;
    const ccomplex_t IT_0122 = m_t*IT_0070;
    const ccomplex_t IT_0123 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0124 = IT_0036*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0127 = IT_0004*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = m_t*conj(N_u4)*e_em*IT_0012*U_st_10;
    const ccomplex_t IT_0130 = IT_0040*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0125 + 3*IT_0128 + 3
      *IT_0131);
    const ccomplex_t IT_0133 = 0.166666666666667*IT_0132;
    const ccomplex_t IT_0134 = IT_0061*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0052*IT_0135;
    const ccomplex_t IT_0137 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0138 = IT_0036*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0141 = IT_0004*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = m_t*conj(N_u4)*e_em*IT_0012*U_st_11;
    const ccomplex_t IT_0144 = IT_0040*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*(IT_0139 + 3*IT_0142 + 3
      *IT_0145);
    const ccomplex_t IT_0147 = 0.166666666666667*IT_0146;
    const ccomplex_t IT_0148 = IT_0047*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0034*IT_0149;
    const ccomplex_t IT_0151 = (-2)*IT_0116 + -IT_0136 + -IT_0150;
    const ccomplex_t IT_0152 = m_t*m_N_4;
    const ccomplex_t IT_0153 = IT_0070*IT_0152;
    const ccomplex_t IT_0154 = -IT_0019 + -IT_0032;
    const ccomplex_t IT_0155 = pow(m_W, 2);
    const ccomplex_t IT_0156 = IT_0152*IT_0155;
    const ccomplex_t IT_0157 = 3*IT_0156;
    const ccomplex_t IT_0158 = m_b*conj(N_d4)*e_em*IT_0012*U_sb_10;
    const ccomplex_t IT_0159 = IT_0075*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0162 = IT_0004*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0165 = IT_0036*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0160 + -IT_0163 +
       0.333333333333333*IT_0166);
    const ccomplex_t IT_0168 = 3*IT_0167;
    const ccomplex_t IT_0169 = 0.166666666666667*IT_0168;
    const ccomplex_t IT_0170 = IT_0089*IT_0169;
    const ccomplex_t IT_0171 = IT_0086*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = s_14*IT_0155;
    const ccomplex_t IT_0174 = IT_0066*IT_0173;
    const ccomplex_t IT_0175 = 0.5*IT_0174;
    const ccomplex_t IT_0176 = s_14 + IT_0175;
    const ccomplex_t IT_0177 = m_t*IT_0176;
    const ccomplex_t IT_0178 = 2*IT_0177;
    const ccomplex_t IT_0179 = pow(m_W, 4);
    const ccomplex_t IT_0180 = s_13*IT_0179;
    const ccomplex_t IT_0181 = IT_0066*IT_0180;
    const ccomplex_t IT_0182 = s_14*s_34;
    const ccomplex_t IT_0183 = 24*IT_0182;
    const ccomplex_t IT_0184 = s_13*IT_0155;
    const ccomplex_t IT_0185 = (-12)*IT_0184;
    const ccomplex_t IT_0186 = 6*IT_0181 + IT_0183 + IT_0185;
    const ccomplex_t IT_0187 = (-3)*IT_0156;
    const ccomplex_t IT_0188 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0189 = IT_0004*IT_0188;
    const ccomplex_t IT_0190 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0191 = IT_0004*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*(IT_0189 + 0.5*IT_0192);
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = IT_0082*IT_0194;
    const ccomplex_t IT_0196 = IT_0002*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = -conj(IT_0116);
    const ccomplex_t IT_0200 = conj(IT_0198) + IT_0199;
    const ccomplex_t IT_0201 = s_14*s_34*IT_0066;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = s_13 + IT_0202;
    const ccomplex_t IT_0204 = s_13*IT_0203;
    const ccomplex_t IT_0205 = -IT_0050 + -IT_0064;
    const ccomplex_t IT_0206 = (-6)*IT_0065;
    const ccomplex_t IT_0207 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0208 = IT_0152*IT_0203;
    const ccomplex_t IT_0209 = pow(s_34, 2);
    const ccomplex_t IT_0210 = IT_0066*IT_0209;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0000 + IT_0211;
    const ccomplex_t IT_0213 = s_13*IT_0212;
    const ccomplex_t IT_0214 = m_t*IT_0203;
    const ccomplex_t IT_0215 = (-6)*IT_0205;
    const ccomplex_t IT_0216 = m_N_4*IT_0203;
    const ccomplex_t IT_0217 = s_34*IT_0155;
    const ccomplex_t IT_0218 = IT_0066*IT_0217;
    const ccomplex_t IT_0219 = (-0.25)*IT_0218;
    const ccomplex_t IT_0220 = s_34 + IT_0219;
    const ccomplex_t IT_0221 = m_N_4*IT_0220;
    const ccomplex_t IT_0222 = (-4)*IT_0221;
    const ccomplex_t IT_0223 = (-6)*IT_0121;
    const ccomplex_t IT_0224 = IT_0136 + IT_0150 + 2*IT_0172;
    const ccomplex_t IT_0225 = (-6)*IT_0151;
    const ccomplex_t IT_0226 = IT_0152*IT_0212;
    const ccomplex_t IT_0227 = (-6)*IT_0224;
    const ccomplex_t IT_0228 = m_N_4*IT_0212;
    const ccomplex_t IT_0229 = (-6)*IT_0181;
    const ccomplex_t IT_0230 = (-24)*IT_0182;
    const ccomplex_t IT_0231 = 12*IT_0184;
    const ccomplex_t IT_0232 = IT_0229 + IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = IT_0183 + IT_0229;
    const ccomplex_t IT_0234 = (-6)*conj(IT_0121);
    const ccomplex_t IT_0235 = conj(IT_0154)*IT_0187;
    const ccomplex_t IT_0236 = m_t*IT_0172;
    const ccomplex_t IT_0237 = m_N_4*IT_0032;
    const ccomplex_t IT_0238 = m_N_4*IT_0019;
    const ccomplex_t IT_0239 = -IT_0236 + -IT_0237 + -IT_0238;
    const ccomplex_t IT_0240 = IT_0029*IT_0106;
    const ccomplex_t IT_0241 = IT_0020*IT_0118*IT_0240;
    const ccomplex_t IT_0242 = IT_0016*IT_0194;
    const ccomplex_t IT_0243 = IT_0002*IT_0072*IT_0242;
    const ccomplex_t IT_0244 = IT_0241 + IT_0243;
    const ccomplex_t IT_0245 = m_N_4*IT_0197;
    const ccomplex_t IT_0246 = 2*IT_0197;
    const ccomplex_t IT_0247 = (-6)*IT_0246;
    const ccomplex_t IT_0248 = 0.166666666666667*IT_0247;
    const ccomplex_t IT_0249 = IT_0151 + IT_0248;
    const ccomplex_t IT_0250 = (-2)*IT_0182;
    const ccomplex_t IT_0251 = IT_0184 + IT_0250;
    const ccomplex_t IT_0252 = IT_0066*IT_0251;
    const ccomplex_t IT_0253 = (-6)*IT_0252;
    const ccomplex_t IT_0254 = 12*s_13;
    const ccomplex_t IT_0255 = IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = IT_0239 + IT_0244;
    const ccomplex_t IT_0257 = 6*IT_0154;
    const ccomplex_t IT_0258 = (-0.25)*IT_0174;
    const ccomplex_t IT_0259 = s_14 + IT_0258;
    const ccomplex_t IT_0260 = m_t*IT_0259;
    const ccomplex_t IT_0261 = (-4)*IT_0260;
    const ccomplex_t IT_0262 = 6*IT_0172;
    const ccomplex_t IT_0263 = -IT_0116;
    const ccomplex_t IT_0264 = IT_0198 + IT_0263;
    const ccomplex_t IT_0265 = (-6)*IT_0264;
    const ccomplex_t IT_0266 = IT_0206*IT_0216 + IT_0214*IT_0227 + IT_0215
      *IT_0228 + IT_0152*((-18)*IT_0121 + (-18)*IT_0245) + 6*IT_0122*IT_0249 +
       IT_0255*IT_0256 + IT_0222*IT_0257 + IT_0261*IT_0262 + IT_0178*IT_0265;
    const ccomplex_t IT_0267 = 0.166666666666667*IT_0225;
    const ccomplex_t IT_0268 = IT_0246 + IT_0267;
    const ccomplex_t IT_0269 = 0.5*IT_0218;
    const ccomplex_t IT_0270 = s_34 + IT_0269;
    const ccomplex_t IT_0271 = m_N_4*IT_0270;
    const ccomplex_t IT_0272 = 2*IT_0271;
    const ccomplex_t IT_0273 = 6*IT_0065*IT_0122 + IT_0214*IT_0215 + IT_0227
      *IT_0228 + IT_0152*((-18)*IT_0239 + (-18)*IT_0244) + (IT_0121 + IT_0245)
      *IT_0255 + IT_0178*IT_0257 + IT_0222*IT_0265 + 6*IT_0216*IT_0268 + IT_0262
      *IT_0272;
    const ccomplex_t IT_0274 = 0.166666666666667*IT_0153;
    const ccomplex_t IT_0275 = (-6)*conj(IT_0151);
    const ccomplex_t IT_0276 = 0.166666666666667*IT_0275;
    const ccomplex_t IT_0277 = conj(IT_0246) + IT_0276;
    const ccomplex_t IT_0278 = (-6)*conj(IT_0205);
    const ccomplex_t IT_0279 = 0.166666666666667*IT_0208;
    const ccomplex_t IT_0280 = 0.166666666666667*IT_0216;
    const ccomplex_t IT_0281 = s_13*s_14;
    const ccomplex_t IT_0282 = s_34*IT_0046;
    const ccomplex_t IT_0283 = -IT_0282;
    const ccomplex_t IT_0284 = IT_0281 + IT_0283;
    const ccomplex_t IT_0285 = (-6)*conj(IT_0246);
    const ccomplex_t IT_0286 = 0.166666666666667*IT_0285;
    const ccomplex_t IT_0287 = conj(IT_0151) + IT_0286;
    const ccomplex_t IT_0288 = IT_0172 + IT_0198;
    const ccomplex_t IT_0289 = 0.166666666666667*IT_0288;
    const ccomplex_t IT_0290 = s_14*IT_0000;
    const ccomplex_t IT_0291 = s_13*s_34;
    const ccomplex_t IT_0292 = -IT_0291;
    const ccomplex_t IT_0293 = IT_0290 + IT_0292;
    const ccomplex_t IT_0294 = (-6)*conj(IT_0224);
    const ccomplex_t IT_0295 = 6*IT_0181 + IT_0230;
    const ccomplex_t IT_0296 = 0.166666666666667*conj(IT_0172);
    const ccomplex_t IT_0297 = 0.166666666666667*conj(IT_0116);
    const ccomplex_t IT_0298 = (-18)*conj(IT_0121);
    const ccomplex_t IT_0299 = -IT_0200;
    const ccomplex_t IT_0300 = 0.166666666666667*IT_0214;
    const ccomplex_t IT_0301 = 0.166666666666667*IT_0228;
    const ccomplex_t IT_0302 = conj(IT_0065)*IT_0122;
    const ccomplex_t IT_0303 = conj(IT_0154)*IT_0178;
    const ccomplex_t IT_0304 = 6*IT_0065*(conj(IT_0065)*IT_0071 + conj(IT_0121
      )*IT_0122 + conj(IT_0151)*IT_0153) + 6*IT_0154*(IT_0157*conj(IT_0172) +
       conj(IT_0121)*IT_0178 + 0.166666666666667*conj(IT_0154)*IT_0186 + 
      -IT_0187*IT_0200) + IT_0204*(conj(IT_0205)*IT_0206 + IT_0205*IT_0207) + 
      (conj(IT_0151)*IT_0208 + conj(IT_0205)*IT_0213 + conj(IT_0121)*IT_0214)
      *IT_0215 + (conj(IT_0151)*IT_0216 + conj(IT_0198)*IT_0222)*IT_0223 + conj
      (IT_0224)*(IT_0206*IT_0208 + IT_0204*IT_0225 + IT_0215*IT_0226 + IT_0213
      *IT_0227 + IT_0223*IT_0228) + IT_0198*(IT_0186*conj(IT_0198) + conj
      (IT_0116)*IT_0232 + conj(IT_0172)*IT_0233 + IT_0222*IT_0234 + (-6)*IT_0235
      ) + (conj(IT_0239) + conj(IT_0244))*IT_0266 + conj(IT_0245)*IT_0273 + 6
      *conj(IT_0246)*(IT_0205*IT_0208 + 0.166666666666667*IT_0071*IT_0225 +
       IT_0206*IT_0274) + 6*IT_0246*(conj(IT_0205)*IT_0208 + conj(IT_0121)
      *IT_0216 + IT_0204*conj(IT_0224) + IT_0207*IT_0274 + IT_0071*IT_0277) + 6
      *IT_0151*(IT_0071*conj(IT_0151) + conj(IT_0065)*IT_0153 + IT_0278*IT_0279 
      + IT_0234*IT_0280) + 6*IT_0284*(conj(IT_0065)*IT_0154 + IT_0065*conj
      (IT_0154) + conj(IT_0116)*IT_0249 + (conj(IT_0172) + conj(IT_0198))
      *IT_0268 + IT_0116*IT_0287 + conj(IT_0246)*IT_0288 + IT_0275*IT_0289) + 6
      *IT_0293*(conj(IT_0154)*IT_0205 + IT_0154*conj(IT_0205) + conj(IT_0116)
      *IT_0224 + IT_0116*conj(IT_0224) + 0.166666666666667*(conj(IT_0172) + conj
      (IT_0198))*IT_0227 + IT_0289*IT_0294) + 6*IT_0116*(conj(IT_0121)*IT_0222 +
       0.166666666666667*conj(IT_0198)*IT_0232 + IT_0235 + IT_0295*IT_0296 +
       IT_0186*IT_0297) + IT_0262*(conj(IT_0154)*IT_0157 + 0.166666666666667
      *conj(IT_0198)*IT_0233 + conj(IT_0121)*IT_0272 + IT_0186*IT_0296 + IT_0295
      *IT_0297) + 6*IT_0224*(IT_0204*IT_0277 + 0.166666666666667*IT_0226*IT_0278
       + IT_0207*IT_0279 + IT_0234*IT_0301) + 6*IT_0256*(conj(IT_0154)*IT_0222 +
       conj(IT_0172)*IT_0261 + IT_0207*IT_0280 + IT_0122*IT_0287 +
       0.166666666666667*IT_0152*IT_0298 + IT_0178*IT_0299 + IT_0294*IT_0300 +
       IT_0278*IT_0301) + 6*IT_0121*(conj(IT_0116)*IT_0222 + IT_0216*conj
      (IT_0246) + 0.166666666666667*conj(IT_0121)*IT_0255 + conj(IT_0172)
      *IT_0272 + IT_0278*IT_0300 + IT_0302 + IT_0303) + 6*IT_0245*
      (0.166666666666667*conj(IT_0121)*IT_0255 + conj(IT_0172)*IT_0272 + IT_0216
      *IT_0277 + IT_0222*IT_0299 + IT_0278*IT_0300 + IT_0294*IT_0301 + IT_0302 +
       IT_0303);
    return create_ccomplex_return(IT_0304);
}

