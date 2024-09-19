#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_sb_1_to_t_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_sb_1_to_t_anti_W(
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
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
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
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0003*IT_0005*IT_0007*U_sb_00;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + (-0.5)*IT_0014);
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = IT_0002*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0021 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0022 = IT_0007*IT_0021;
    const ccomplex_t IT_0023 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0024 = IT_0007*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0022 + (-0.5)*IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0003*IT_0005*IT_0007*U_sb_00;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0007*U_sb_00*conj(U_st_00);
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = cos(theta_W);
    const ccomplex_t IT_0036 = cpow(IT_0035, -1);
    const ccomplex_t IT_0037 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = IT_0005*IT_0007;
    const ccomplex_t IT_0041 = m_t*N_u1*e_em*IT_0003*U_st_00;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0039 + (-0.75)*IT_0043);
    const ccomplex_t IT_0045 = 0.666666666666667*IT_0044;
    const ccomplex_t IT_0046 = pow(m_N_1, 2);
    const ccomplex_t IT_0047 = cpow((-2)*s_13 + IT_0000 + IT_0046 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0034*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0007*U_sb_00*conj(U_st_01);
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0054 = IT_0036*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = m_t*N_u1*e_em*IT_0003*U_st_01;
    const ccomplex_t IT_0057 = IT_0040*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + (-0.75)*IT_0058);
    const ccomplex_t IT_0060 = 0.666666666666667*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_13 + IT_0000 + IT_0046 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
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
    const ccomplex_t IT_0072 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0073 = IT_0036*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0076 = IT_0007*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = m_t*conj(N_u1)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0079 = IT_0040*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0074 + 3*IT_0077 + 3
      *IT_0080);
    const ccomplex_t IT_0082 = 0.166666666666667*IT_0081;
    const ccomplex_t IT_0083 = IT_0047*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = IT_0034*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0087 = IT_0036*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0090 = IT_0007*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = m_t*conj(N_u1)*e_em*IT_0003*U_st_11;
    const ccomplex_t IT_0093 = IT_0040*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0088 + 3*IT_0091 + 3
      *IT_0094);
    const ccomplex_t IT_0096 = 0.166666666666667*IT_0095;
    const ccomplex_t IT_0097 = IT_0061*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0052*IT_0098;
    const ccomplex_t IT_0100 = -IT_0085 + -IT_0099;
    const ccomplex_t IT_0101 = m_t*m_N_1;
    const ccomplex_t IT_0102 = IT_0070*IT_0101;
    const ccomplex_t IT_0103 = -IT_0019 + -IT_0032;
    const ccomplex_t IT_0104 = s_13*s_14;
    const ccomplex_t IT_0105 = s_34*IT_0046;
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = IT_0104 + IT_0106;
    const ccomplex_t IT_0108 = pow(m_b, 2);
    const ccomplex_t IT_0109 = cpow((-2)*s_12 + -IT_0001 + -IT_0046 + IT_0108 
      + -reg_prop, -1);
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0007;
    const ccomplex_t IT_0111 = 0.5*IT_0110;
    const ccomplex_t IT_0112 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0113 = IT_0036*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0116 = IT_0007*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = cos(beta);
    const ccomplex_t IT_0119 = cpow(IT_0118, -1);
    const ccomplex_t IT_0120 = IT_0007*IT_0119;
    const ccomplex_t IT_0121 = m_b*conj(N_d1)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0114 + (-3)*IT_0117 + 3
      *IT_0123);
    const ccomplex_t IT_0125 = (-0.333333333333333)*IT_0124;
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = 3*IT_0126;
    const ccomplex_t IT_0128 = 0.166666666666667*IT_0127;
    const ccomplex_t IT_0129 = IT_0111*IT_0128;
    const ccomplex_t IT_0130 = IT_0109*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = pow(m_W, 2);
    const ccomplex_t IT_0133 = IT_0101*IT_0132;
    const ccomplex_t IT_0134 = 3*IT_0133;
    const ccomplex_t IT_0135 = m_b*conj(U_d2)*V_tb*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0136 = IT_0120*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0139 = IT_0007*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0137 + (-2)*IT_0139);
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0143 = IT_0007*IT_0142;
    const ccomplex_t IT_0144 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0145 = IT_0007*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0143 + 0.5*IT_0146);
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = IT_0141*IT_0148;
    const ccomplex_t IT_0150 = IT_0020*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = m_b*conj(U_d1)*V_tb*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0153 = IT_0120*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0156 = IT_0007*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0154 + (-2)*IT_0156);
    const ccomplex_t IT_0158 = (-0.5)*IT_0157;
    const ccomplex_t IT_0159 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0160 = IT_0007*IT_0159;
    const ccomplex_t IT_0161 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0162 = IT_0007*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0160 + 0.5*IT_0163);
    const ccomplex_t IT_0165 = -IT_0164;
    const ccomplex_t IT_0166 = IT_0158*IT_0165;
    const ccomplex_t IT_0167 = IT_0002*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = -IT_0151 + -IT_0168;
    const ccomplex_t IT_0170 = (-3)*IT_0133;
    const ccomplex_t IT_0171 = -IT_0050 + -IT_0064;
    const ccomplex_t IT_0172 = s_14*IT_0000;
    const ccomplex_t IT_0173 = s_13*s_34;
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = IT_0172 + IT_0174;
    const ccomplex_t IT_0176 = pow(m_W, 4);
    const ccomplex_t IT_0177 = s_13*IT_0176;
    const ccomplex_t IT_0178 = IT_0066*IT_0177;
    const ccomplex_t IT_0179 = s_14*s_34;
    const ccomplex_t IT_0180 = 24*IT_0179;
    const ccomplex_t IT_0181 = s_13*IT_0132;
    const ccomplex_t IT_0182 = (-12)*IT_0181;
    const ccomplex_t IT_0183 = 6*IT_0178 + IT_0180 + IT_0182;
    const ccomplex_t IT_0184 = s_14*IT_0132;
    const ccomplex_t IT_0185 = IT_0066*IT_0184;
    const ccomplex_t IT_0186 = 0.5*IT_0185;
    const ccomplex_t IT_0187 = s_14 + IT_0186;
    const ccomplex_t IT_0188 = m_t*IT_0187;
    const ccomplex_t IT_0189 = 2*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0191 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0192 = IT_0036*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = m_b*N_d1*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0195 = IT_0120*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + 1.5*IT_0196);
    const ccomplex_t IT_0198 = (-0.333333333333333)*IT_0197;
    const ccomplex_t IT_0199 = IT_0111*IT_0198;
    const ccomplex_t IT_0200 = IT_0109*IT_0190*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0202 = IT_0016*IT_0158;
    const ccomplex_t IT_0203 = IT_0002*IT_0201*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0205 = IT_0027*IT_0141;
    const ccomplex_t IT_0206 = IT_0020*IT_0204*IT_0205;
    const ccomplex_t IT_0207 = -IT_0200 + -IT_0203 + -IT_0206;
    const ccomplex_t IT_0208 = IT_0151 + IT_0168;
    const ccomplex_t IT_0209 = m_N_1*IT_0208;
    const ccomplex_t IT_0210 = IT_0207 + IT_0209;
    const ccomplex_t IT_0211 = s_34*IT_0132;
    const ccomplex_t IT_0212 = IT_0066*IT_0211;
    const ccomplex_t IT_0213 = (-0.25)*IT_0212;
    const ccomplex_t IT_0214 = s_34 + IT_0213;
    const ccomplex_t IT_0215 = m_N_1*IT_0214;
    const ccomplex_t IT_0216 = (-4)*IT_0215;
    const ccomplex_t IT_0217 = m_N_1*IT_0019;
    const ccomplex_t IT_0218 = m_N_1*IT_0032;
    const ccomplex_t IT_0219 = m_t*IT_0131;
    const ccomplex_t IT_0220 = -IT_0217 + -IT_0218 + -IT_0219;
    const ccomplex_t IT_0221 = IT_0029*IT_0148;
    const ccomplex_t IT_0222 = IT_0020*IT_0204*IT_0221;
    const ccomplex_t IT_0223 = IT_0009*IT_0165;
    const ccomplex_t IT_0224 = IT_0002*IT_0201*IT_0223;
    const ccomplex_t IT_0225 = IT_0222 + IT_0224;
    const ccomplex_t IT_0226 = IT_0220 + IT_0225;
    const ccomplex_t IT_0227 = (-18)*conj(IT_0207);
    const ccomplex_t IT_0228 = s_14*s_34*IT_0066;
    const ccomplex_t IT_0229 = -IT_0228;
    const ccomplex_t IT_0230 = s_13 + IT_0229;
    const ccomplex_t IT_0231 = s_13*IT_0230;
    const ccomplex_t IT_0232 = 2*IT_0151 + 2*IT_0168;
    const ccomplex_t IT_0233 = IT_0085 + IT_0099 + 2*IT_0131;
    const ccomplex_t IT_0234 = IT_0101*IT_0230;
    const ccomplex_t IT_0235 = m_N_1*IT_0230;
    const ccomplex_t IT_0236 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0237 = 0.5*IT_0212;
    const ccomplex_t IT_0238 = s_34 + IT_0237;
    const ccomplex_t IT_0239 = m_N_1*IT_0238;
    const ccomplex_t IT_0240 = 2*IT_0239;
    const ccomplex_t IT_0241 = (-2)*IT_0179;
    const ccomplex_t IT_0242 = IT_0181 + IT_0241;
    const ccomplex_t IT_0243 = IT_0066*IT_0242;
    const ccomplex_t IT_0244 = (-6)*IT_0243;
    const ccomplex_t IT_0245 = 12*s_13;
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = m_t*IT_0230;
    const ccomplex_t IT_0248 = (-6)*IT_0171;
    const ccomplex_t IT_0249 = pow(s_34, 2);
    const ccomplex_t IT_0250 = IT_0066*IT_0249;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = IT_0000 + IT_0251;
    const ccomplex_t IT_0253 = IT_0101*IT_0252;
    const ccomplex_t IT_0254 = m_N_1*IT_0252;
    const ccomplex_t IT_0255 = (-6)*conj(IT_0171);
    const ccomplex_t IT_0256 = (-6)*IT_0178;
    const ccomplex_t IT_0257 = IT_0180 + IT_0256;
    const ccomplex_t IT_0258 = (-6)*IT_0100;
    const ccomplex_t IT_0259 = 0.166666666666667*IT_0258;
    const ccomplex_t IT_0260 = IT_0232 + IT_0259;
    const ccomplex_t IT_0261 = (-0.25)*IT_0185;
    const ccomplex_t IT_0262 = s_14 + IT_0261;
    const ccomplex_t IT_0263 = m_t*IT_0262;
    const ccomplex_t IT_0264 = (-4)*IT_0263;
    const ccomplex_t IT_0265 = (-6)*IT_0233;
    const ccomplex_t IT_0266 = IT_0131 + IT_0169;
    const ccomplex_t IT_0267 = (-6)*conj(IT_0233);
    const ccomplex_t IT_0268 = (-6)*IT_0207;
    const ccomplex_t IT_0269 = (-6)*conj(IT_0169)*(IT_0103*IT_0170 + (
      -0.166666666666667)*IT_0169*IT_0183 + IT_0189*IT_0226 + (
      -0.166666666666667)*IT_0131*IT_0257 + -IT_0107*IT_0260 + (
      -0.166666666666667)*IT_0175*IT_0265 + IT_0216*(IT_0209 + (
      -0.166666666666667)*IT_0268));
    const ccomplex_t IT_0270 = (-6)*conj(IT_0207);
    const ccomplex_t IT_0271 = (-6)*conj(IT_0232);
    const ccomplex_t IT_0272 = m_t*IT_0070;
    const ccomplex_t IT_0273 = 6*IT_0272;
    const ccomplex_t IT_0274 = (-6)*IT_0065;
    const ccomplex_t IT_0275 = 6*IT_0103;
    const ccomplex_t IT_0276 = 6*IT_0131;
    const ccomplex_t IT_0277 = (-6)*IT_0232;
    const ccomplex_t IT_0278 = (-6)*IT_0169;
    const ccomplex_t IT_0279 = IT_0101*((-18)*IT_0207 + (-18)*IT_0209) +
       IT_0226*IT_0246 + IT_0248*IT_0254 + IT_0247*IT_0265 + IT_0235*IT_0274 +
       IT_0216*IT_0275 + IT_0264*IT_0276 + IT_0273*(IT_0100 + 0.166666666666667
      *IT_0277) + IT_0189*IT_0278;
    const ccomplex_t IT_0280 = IT_0101*((-18)*IT_0220 + (-18)*IT_0225) +
       IT_0210*IT_0246 + IT_0247*IT_0248 + 6*IT_0235*IT_0260 + IT_0254*IT_0265 +
       6*IT_0065*IT_0272 + IT_0189*IT_0275 + IT_0240*IT_0276 + IT_0216*IT_0278;
    const ccomplex_t IT_0281 = 0.166666666666667*IT_0234;
    const ccomplex_t IT_0282 = s_13*IT_0252;
    const ccomplex_t IT_0283 = 0.166666666666667*IT_0282;
    const ccomplex_t IT_0284 = (-6)*conj(IT_0100);
    const ccomplex_t IT_0285 = IT_0107*IT_0266;
    const ccomplex_t IT_0286 = 6*conj(IT_0065)*(IT_0065*IT_0071 + IT_0100
      *IT_0102 + IT_0103*IT_0107) + 6*conj(IT_0103)*(IT_0065*IT_0107 + IT_0131
      *IT_0134 + -IT_0169*IT_0170 + IT_0171*IT_0175 + 0.166666666666667*IT_0103
      *IT_0183 + IT_0189*IT_0210 + IT_0216*IT_0226) + IT_0101*IT_0226*IT_0227 + 
      (IT_0171*IT_0231 + IT_0102*IT_0232 + IT_0233*IT_0234 + IT_0226*IT_0235)
      *IT_0236 + 6*conj(IT_0207)*(IT_0103*IT_0189 + IT_0232*IT_0235 + IT_0131
      *IT_0240 + 0.166666666666667*IT_0210*IT_0246 + 0.166666666666667*IT_0247
      *IT_0248) + (IT_0100*IT_0234 + IT_0210*IT_0247 + IT_0233*IT_0253 + IT_0226
      *IT_0254)*IT_0255 + 6*conj(IT_0131)*(IT_0103*IT_0134 + 0.166666666666667
      *IT_0131*IT_0183 + IT_0210*IT_0240 + 0.166666666666667*IT_0169*IT_0257 +
       IT_0107*IT_0260 + IT_0226*IT_0264 + 0.166666666666667*IT_0175*IT_0265) + 
      (IT_0226*IT_0247 + IT_0209*IT_0254 + IT_0175*IT_0266)*IT_0267 + IT_0269 + 
      (IT_0169*IT_0216 + IT_0100*IT_0235 + IT_0233*IT_0254)*IT_0270 + (IT_0065
      *conj(IT_0207) + conj(IT_0065)*IT_0210 + IT_0226*(conj(IT_0100) +
       0.166666666666667*IT_0271))*IT_0273 + (conj(IT_0220) + conj(IT_0225))
      *IT_0279 + conj(IT_0209)*IT_0280 + 6*conj(IT_0100)*(IT_0071*IT_0100 +
       IT_0065*IT_0102 + 0.166666666666667*IT_0235*IT_0268 + IT_0248*IT_0281) +
       6*conj(IT_0171)*(IT_0103*IT_0175 + IT_0232*IT_0234 + 0.166666666666667
      *IT_0231*IT_0274 + IT_0248*IT_0283) + 6*conj(IT_0233)*(0.166666666666667
      *IT_0248*IT_0253 + IT_0231*IT_0260 + 0.166666666666667*IT_0254*IT_0268 +
       IT_0274*IT_0281 + IT_0265*IT_0283) + IT_0284*(IT_0071*IT_0232 + IT_0231
      *IT_0233 + IT_0209*IT_0235 + IT_0285) + 6*conj(IT_0232)*(IT_0231*IT_0233 +
       IT_0171*IT_0234 + IT_0210*IT_0235 + IT_0071*IT_0260 + 0.166666666666667
      *IT_0102*IT_0274 + IT_0285);
    return create_ccomplex_return(IT_0286);
}

