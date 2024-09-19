#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_st_1_to_b_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_st_1_to_b_W(
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
    const creal_t m_st_1 = param->m_st_1;
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
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_01)*U_st_00;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0001*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = m_b*N_d2*e_em*IT_0012*U_sb_01;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0008 + 1.5*IT_0015);
    const ccomplex_t IT_0017 = (-0.333333333333333)*IT_0016;
    const ccomplex_t IT_0018 = pow(m_b, 2);
    const ccomplex_t IT_0019 = pow(m_N_2, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + 0.5*IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0001*IT_0010*IT_0012*U_st_00;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = pow(m_st_1, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0018 + IT_0034 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0039 = IT_0001*IT_0038;
    const ccomplex_t IT_0040 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + 0.5*IT_0042);
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0001*IT_0010*IT_0012*U_st_00;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0018 + IT_0034 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0054 = IT_0005*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = m_b*N_d2*e_em*IT_0012*U_sb_00;
    const ccomplex_t IT_0057 = IT_0011*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + 1.5*IT_0058);
    const ccomplex_t IT_0060 = (-0.333333333333333)*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0063;
    const ccomplex_t IT_0065 = IT_0023 + (-2)*IT_0037 + (-2)*IT_0050 + IT_0064;
    const ccomplex_t IT_0066 = pow(m_W, -2);
    const ccomplex_t IT_0067 = pow(s_14, 2);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = -IT_0019;
    const ccomplex_t IT_0070 = IT_0068 + IT_0069;
    const ccomplex_t IT_0071 = s_13*IT_0070;
    const ccomplex_t IT_0072 = m_b*conj(N_d2)*e_em*IT_0012*U_sb_11;
    const ccomplex_t IT_0073 = IT_0011*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0076 = IT_0001*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0079 = IT_0005*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0074 + -IT_0077 +
       0.333333333333333*IT_0080);
    const ccomplex_t IT_0082 = 3*IT_0081;
    const ccomplex_t IT_0083 = 0.166666666666667*IT_0082;
    const ccomplex_t IT_0084 = IT_0020*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0003*IT_0085;
    const ccomplex_t IT_0087 = m_b*conj(N_d2)*e_em*IT_0012*U_sb_10;
    const ccomplex_t IT_0088 = IT_0011*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0091 = IT_0005*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0094 = IT_0001*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0089 + 0.333333333333333
      *IT_0092 + -IT_0095);
    const ccomplex_t IT_0097 = 3*IT_0096;
    const ccomplex_t IT_0098 = 0.166666666666667*IT_0097;
    const ccomplex_t IT_0099 = IT_0061*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0052*IT_0100;
    const ccomplex_t IT_0102 = -IT_0086 + -IT_0101;
    const ccomplex_t IT_0103 = m_b*m_N_2;
    const ccomplex_t IT_0104 = IT_0070*IT_0103;
    const ccomplex_t IT_0105 = IT_0037 + IT_0050;
    const ccomplex_t IT_0106 = s_13*s_14;
    const ccomplex_t IT_0107 = s_34*IT_0019;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = IT_0106 + IT_0108;
    const ccomplex_t IT_0110 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0111 = IT_0001*IT_0110;
    const ccomplex_t IT_0112 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0113 = IT_0001*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + (-0.5)*IT_0114);
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = V_tb*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0118 = IT_0001*IT_0117;
    const ccomplex_t IT_0119 = sin(beta);
    const ccomplex_t IT_0120 = cpow(IT_0119, -1);
    const ccomplex_t IT_0121 = IT_0001*IT_0120;
    const ccomplex_t IT_0122 = m_t*V_tb*conj(V_u1)*e_em*IT_0012*U_st_10;
    const ccomplex_t IT_0123 = IT_0121*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0118 + (-0.5)*IT_0124);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = IT_0116*IT_0126;
    const ccomplex_t IT_0128 = IT_0035*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0131 = IT_0001*IT_0130;
    const ccomplex_t IT_0132 = m_t*V_tb*conj(V_u2)*e_em*IT_0012*U_st_10;
    const ccomplex_t IT_0133 = IT_0121*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0131 + (-0.5)*IT_0134);
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0138 = IT_0001*IT_0137;
    const ccomplex_t IT_0139 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0140 = IT_0001*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0138 + (-0.5)*IT_0141);
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = IT_0136*IT_0143;
    const ccomplex_t IT_0145 = IT_0048*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = -IT_0129 + -IT_0146;
    const ccomplex_t IT_0148 = pow(m_W, 2);
    const ccomplex_t IT_0149 = IT_0103*IT_0148;
    const ccomplex_t IT_0150 = (-3)*IT_0149;
    const ccomplex_t IT_0151 = 3*IT_0149;
    const ccomplex_t IT_0152 = cpow((-2)*s_12 + -IT_0019 + -IT_0034 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0153 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0154 = IT_0005*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0157 = IT_0001*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = m_t*conj(N_u2)*e_em*IT_0012*U_st_10;
    const ccomplex_t IT_0160 = IT_0121*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*(IT_0155 + 3*IT_0158 + 3
      *IT_0161);
    const ccomplex_t IT_0163 = 0.166666666666667*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001;
    const ccomplex_t IT_0165 = 0.5*IT_0164;
    const ccomplex_t IT_0166 = IT_0163*IT_0165;
    const ccomplex_t IT_0167 = IT_0152*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = -IT_0023 + -IT_0064;
    const ccomplex_t IT_0170 = s_14*IT_0018;
    const ccomplex_t IT_0171 = s_13*s_34;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = IT_0170 + IT_0172;
    const ccomplex_t IT_0174 = pow(m_W, 4);
    const ccomplex_t IT_0175 = s_13*IT_0174;
    const ccomplex_t IT_0176 = IT_0066*IT_0175;
    const ccomplex_t IT_0177 = s_14*s_34;
    const ccomplex_t IT_0178 = 24*IT_0177;
    const ccomplex_t IT_0179 = s_13*IT_0148;
    const ccomplex_t IT_0180 = (-12)*IT_0179;
    const ccomplex_t IT_0181 = 6*IT_0176 + IT_0178 + IT_0180;
    const ccomplex_t IT_0182 = s_14*IT_0148;
    const ccomplex_t IT_0183 = IT_0066*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0183;
    const ccomplex_t IT_0185 = s_14 + IT_0184;
    const ccomplex_t IT_0186 = m_b*IT_0185;
    const ccomplex_t IT_0187 = 2*IT_0186;
    const ccomplex_t IT_0188 = IT_0030*IT_0126;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0190 = IT_0035*IT_0188*IT_0189;
    const ccomplex_t IT_0191 = IT_0044*IT_0136;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0193 = IT_0048*IT_0191*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0195 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0196 = IT_0005*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = m_t*N_u2*e_em*IT_0012*U_st_00;
    const ccomplex_t IT_0199 = IT_0121*IT_0198;
    const ccomplex_t IT_0200 = 1.4142135623731*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*(IT_0197 + (-0.75)*IT_0200);
    const ccomplex_t IT_0202 = 0.666666666666667*IT_0201;
    const ccomplex_t IT_0203 = IT_0165*IT_0202;
    const ccomplex_t IT_0204 = IT_0152*IT_0194*IT_0203;
    const ccomplex_t IT_0205 = -IT_0190 + -IT_0193 + -IT_0204;
    const ccomplex_t IT_0206 = IT_0129 + IT_0146;
    const ccomplex_t IT_0207 = m_N_2*IT_0206;
    const ccomplex_t IT_0208 = IT_0205 + IT_0207;
    const ccomplex_t IT_0209 = s_34*IT_0148;
    const ccomplex_t IT_0210 = IT_0066*IT_0209;
    const ccomplex_t IT_0211 = (-0.25)*IT_0210;
    const ccomplex_t IT_0212 = s_34 + IT_0211;
    const ccomplex_t IT_0213 = m_N_2*IT_0212;
    const ccomplex_t IT_0214 = (-4)*IT_0213;
    const ccomplex_t IT_0215 = m_N_2*IT_0037;
    const ccomplex_t IT_0216 = m_N_2*IT_0050;
    const ccomplex_t IT_0217 = m_b*IT_0168;
    const ccomplex_t IT_0218 = IT_0215 + IT_0216 + -IT_0217;
    const ccomplex_t IT_0219 = IT_0032*IT_0116;
    const ccomplex_t IT_0220 = IT_0035*IT_0189*IT_0219;
    const ccomplex_t IT_0221 = IT_0046*IT_0143;
    const ccomplex_t IT_0222 = IT_0048*IT_0192*IT_0221;
    const ccomplex_t IT_0223 = -IT_0220 + -IT_0222;
    const ccomplex_t IT_0224 = IT_0218 + IT_0223;
    const ccomplex_t IT_0225 = 6*conj(IT_0105)*(IT_0065*IT_0109 + -IT_0147
      *IT_0150 + IT_0151*IT_0168 + IT_0169*IT_0173 + 0.166666666666667*IT_0105
      *IT_0181 + IT_0187*IT_0208 + IT_0214*IT_0224);
    const ccomplex_t IT_0226 = (-18)*conj(IT_0205);
    const ccomplex_t IT_0227 = 2*IT_0129 + 2*IT_0146;
    const ccomplex_t IT_0228 = s_14*s_34*IT_0066;
    const ccomplex_t IT_0229 = -IT_0228;
    const ccomplex_t IT_0230 = s_13 + IT_0229;
    const ccomplex_t IT_0231 = m_N_2*IT_0230;
    const ccomplex_t IT_0232 = 0.5*IT_0210;
    const ccomplex_t IT_0233 = s_34 + IT_0232;
    const ccomplex_t IT_0234 = m_N_2*IT_0233;
    const ccomplex_t IT_0235 = 2*IT_0234;
    const ccomplex_t IT_0236 = (-2)*IT_0177;
    const ccomplex_t IT_0237 = IT_0179 + IT_0236;
    const ccomplex_t IT_0238 = IT_0066*IT_0237;
    const ccomplex_t IT_0239 = (-6)*IT_0238;
    const ccomplex_t IT_0240 = 12*s_13;
    const ccomplex_t IT_0241 = IT_0239 + IT_0240;
    const ccomplex_t IT_0242 = m_b*IT_0230;
    const ccomplex_t IT_0243 = (-6)*IT_0169;
    const ccomplex_t IT_0244 = IT_0103*IT_0230;
    const ccomplex_t IT_0245 = IT_0086 + IT_0101 + 2*IT_0168;
    const ccomplex_t IT_0246 = pow(s_34, 2);
    const ccomplex_t IT_0247 = IT_0066*IT_0246;
    const ccomplex_t IT_0248 = -IT_0247;
    const ccomplex_t IT_0249 = IT_0018 + IT_0248;
    const ccomplex_t IT_0250 = IT_0103*IT_0249;
    const ccomplex_t IT_0251 = m_N_2*IT_0249;
    const ccomplex_t IT_0252 = (-6)*conj(IT_0169);
    const ccomplex_t IT_0253 = s_13*IT_0230;
    const ccomplex_t IT_0254 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0255 = (-6)*IT_0176;
    const ccomplex_t IT_0256 = IT_0178 + IT_0255;
    const ccomplex_t IT_0257 = (-6)*IT_0102;
    const ccomplex_t IT_0258 = 0.166666666666667*IT_0257;
    const ccomplex_t IT_0259 = IT_0227 + IT_0258;
    const ccomplex_t IT_0260 = (-0.25)*IT_0183;
    const ccomplex_t IT_0261 = s_14 + IT_0260;
    const ccomplex_t IT_0262 = m_b*IT_0261;
    const ccomplex_t IT_0263 = (-4)*IT_0262;
    const ccomplex_t IT_0264 = (-6)*IT_0245;
    const ccomplex_t IT_0265 = IT_0147 + IT_0168;
    const ccomplex_t IT_0266 = (-6)*conj(IT_0245);
    const ccomplex_t IT_0267 = (-6)*IT_0205;
    const ccomplex_t IT_0268 = (-6)*conj(IT_0147)*(IT_0105*IT_0150 + (
      -0.166666666666667)*IT_0147*IT_0181 + IT_0187*IT_0224 + (
      -0.166666666666667)*IT_0168*IT_0256 + -IT_0109*IT_0259 + (
      -0.166666666666667)*IT_0173*IT_0264 + IT_0214*(IT_0207 + (
      -0.166666666666667)*IT_0267));
    const ccomplex_t IT_0269 = (-6)*conj(IT_0205);
    const ccomplex_t IT_0270 = (-6)*conj(IT_0227);
    const ccomplex_t IT_0271 = m_b*IT_0070;
    const ccomplex_t IT_0272 = 6*IT_0271;
    const ccomplex_t IT_0273 = 6*IT_0105;
    const ccomplex_t IT_0274 = 6*IT_0168;
    const ccomplex_t IT_0275 = (-6)*IT_0147*IT_0214 + IT_0103*((-18)*IT_0218 +
       (-18)*IT_0223) + IT_0208*IT_0241 + IT_0242*IT_0243 + 6*IT_0231*IT_0259 +
       IT_0251*IT_0264 + 6*IT_0065*IT_0271 + IT_0187*IT_0273 + IT_0235*IT_0274;
    const ccomplex_t IT_0276 = (-6)*IT_0065;
    const ccomplex_t IT_0277 = (-6)*IT_0227;
    const ccomplex_t IT_0278 = (-6)*IT_0147*IT_0187 + IT_0103*((-18)*IT_0205 +
       (-18)*IT_0207) + IT_0224*IT_0241 + IT_0243*IT_0251 + IT_0242*IT_0264 +
       IT_0214*IT_0273 + IT_0263*IT_0274 + IT_0231*IT_0276 + IT_0272*(IT_0102 +
       0.166666666666667*IT_0277);
    const ccomplex_t IT_0279 = 0.166666666666667*IT_0244;
    const ccomplex_t IT_0280 = s_13*IT_0249;
    const ccomplex_t IT_0281 = 0.166666666666667*IT_0280;
    const ccomplex_t IT_0282 = (-6)*conj(IT_0102);
    const ccomplex_t IT_0283 = IT_0109*IT_0265;
    const ccomplex_t IT_0284 = 6*conj(IT_0065)*(IT_0065*IT_0071 + IT_0102
      *IT_0104 + IT_0105*IT_0109) + IT_0225 + IT_0103*IT_0224*IT_0226 + 6*conj
      (IT_0205)*(IT_0105*IT_0187 + IT_0227*IT_0231 + IT_0168*IT_0235 +
       0.166666666666667*IT_0208*IT_0241 + 0.166666666666667*IT_0242*IT_0243) + 
      (IT_0208*IT_0242 + IT_0102*IT_0244 + IT_0245*IT_0250 + IT_0224*IT_0251)
      *IT_0252 + (IT_0104*IT_0227 + IT_0224*IT_0231 + IT_0244*IT_0245 + IT_0169
      *IT_0253)*IT_0254 + 6*conj(IT_0168)*(IT_0105*IT_0151 + 0.166666666666667
      *IT_0168*IT_0181 + IT_0208*IT_0235 + 0.166666666666667*IT_0147*IT_0256 +
       IT_0109*IT_0259 + IT_0224*IT_0263 + 0.166666666666667*IT_0173*IT_0264) + 
      (IT_0224*IT_0242 + IT_0207*IT_0251 + IT_0173*IT_0265)*IT_0266 + IT_0268 + 
      (IT_0147*IT_0214 + IT_0102*IT_0231 + IT_0245*IT_0251)*IT_0269 + (IT_0065
      *conj(IT_0205) + conj(IT_0065)*IT_0208 + IT_0224*(conj(IT_0102) +
       0.166666666666667*IT_0270))*IT_0272 + conj(IT_0207)*IT_0275 + (conj
      (IT_0218) + conj(IT_0223))*IT_0278 + 6*conj(IT_0102)*(IT_0071*IT_0102 +
       IT_0065*IT_0104 + 0.166666666666667*IT_0231*IT_0267 + IT_0243*IT_0279) +
       6*conj(IT_0169)*(IT_0105*IT_0173 + IT_0227*IT_0244 + 0.166666666666667
      *IT_0253*IT_0276 + IT_0243*IT_0281) + 6*conj(IT_0245)*(0.166666666666667
      *IT_0243*IT_0250 + IT_0253*IT_0259 + 0.166666666666667*IT_0251*IT_0267 +
       IT_0276*IT_0279 + IT_0264*IT_0281) + IT_0282*(IT_0071*IT_0227 + IT_0207
      *IT_0231 + IT_0245*IT_0253 + IT_0283) + 6*conj(IT_0227)*(IT_0208*IT_0231 +
       IT_0169*IT_0244 + IT_0245*IT_0253 + IT_0071*IT_0259 + 0.166666666666667
      *IT_0104*IT_0276 + IT_0283);
    return create_ccomplex_return(IT_0284);
}

