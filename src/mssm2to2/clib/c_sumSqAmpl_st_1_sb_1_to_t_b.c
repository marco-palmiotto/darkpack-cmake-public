#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_sb_1_to_t_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_sb_1_to_t_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
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
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_st_1, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-6)*IT_0001;
    const ccomplex_t IT_0003 = s_13*s_14;
    const ccomplex_t IT_0004 = 12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_t, 2);
    const ccomplex_t IT_0007 = pow(m_sG, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -IT_0007 +
       reg_prop, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = sin(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = pow(m_W, -1);
    const ccomplex_t IT_0026 = m_b*N_d1*e_em*IT_0025*U_sb_00;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0019 + 1.5*IT_0028);
    const ccomplex_t IT_0030 = (-0.333333333333333)*IT_0029;
    const ccomplex_t IT_0031 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0032 = IT_0023*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = sin(beta);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = IT_0023*IT_0035;
    const ccomplex_t IT_0037 = m_t*conj(N_u1)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0041 = IT_0016*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0033 + IT_0039 +
       0.333333333333333*IT_0042);
    const ccomplex_t IT_0044 = 3*IT_0043;
    const ccomplex_t IT_0045 = 0.166666666666667*IT_0044;
    const ccomplex_t IT_0046 = IT_0030*IT_0045;
    const ccomplex_t IT_0047 = IT_0014*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0050 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0051 = IT_0016*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0054 = IT_0023*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = m_t*conj(N_u2)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0057 = IT_0036*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0052 + 3*IT_0055 + 3
      *IT_0058);
    const ccomplex_t IT_0060 = 0.166666666666667*IT_0059;
    const ccomplex_t IT_0061 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0062 = IT_0016*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = m_b*N_d2*e_em*IT_0025*U_sb_00;
    const ccomplex_t IT_0065 = IT_0024*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + 1.5*IT_0066);
    const ccomplex_t IT_0068 = 0.666666666666667*IT_0067;
    const ccomplex_t IT_0069 = 1.5*IT_0068;
    const ccomplex_t IT_0070 = 0.666666666666667*IT_0069;
    const ccomplex_t IT_0071 = 1.5*IT_0070;
    const ccomplex_t IT_0072 = 0.666666666666667*IT_0071;
    const ccomplex_t IT_0073 = 1.5*IT_0072;
    const ccomplex_t IT_0074 = 0.666666666666667*IT_0073;
    const ccomplex_t IT_0075 = 1.5*IT_0074;
    const ccomplex_t IT_0076 = 0.666666666666667*IT_0075;
    const ccomplex_t IT_0077 = 1.5*IT_0076;
    const ccomplex_t IT_0078 = 0.666666666666667*IT_0077;
    const ccomplex_t IT_0079 = 1.5*IT_0078;
    const ccomplex_t IT_0080 = 0.666666666666667*IT_0079;
    const ccomplex_t IT_0081 = 1.5*IT_0080;
    const ccomplex_t IT_0082 = (-0.333333333333333)*IT_0081;
    const ccomplex_t IT_0083 = IT_0060*IT_0082;
    const ccomplex_t IT_0084 = IT_0049*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0087 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0088 = IT_0023*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = m_t*conj(N_u3)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0091 = IT_0036*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0094 = IT_0016*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0089 + IT_0092 +
       0.333333333333333*IT_0095);
    const ccomplex_t IT_0097 = 3*IT_0096;
    const ccomplex_t IT_0098 = 0.166666666666667*IT_0097;
    const ccomplex_t IT_0099 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0100 = IT_0016*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = m_b*N_d3*e_em*IT_0025*U_sb_00;
    const ccomplex_t IT_0103 = IT_0024*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + 1.5*IT_0104);
    const ccomplex_t IT_0106 = 0.666666666666667*IT_0105;
    const ccomplex_t IT_0107 = 1.5*IT_0106;
    const ccomplex_t IT_0108 = (-0.333333333333333)*IT_0107;
    const ccomplex_t IT_0109 = IT_0098*IT_0108;
    const ccomplex_t IT_0110 = IT_0086*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0113 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0114 = IT_0016*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0117 = IT_0023*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = m_t*conj(N_u4)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0120 = IT_0036*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0115 + 3*IT_0118 + 3
      *IT_0121);
    const ccomplex_t IT_0123 = 0.166666666666667*IT_0122;
    const ccomplex_t IT_0124 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0125 = IT_0016*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = m_b*N_d4*e_em*IT_0025*U_sb_00;
    const ccomplex_t IT_0128 = IT_0024*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + 1.5*IT_0129);
    const ccomplex_t IT_0131 = 0.666666666666667*IT_0130;
    const ccomplex_t IT_0132 = 1.5*IT_0131;
    const ccomplex_t IT_0133 = (-0.333333333333333)*IT_0132;
    const ccomplex_t IT_0134 = IT_0123*IT_0133;
    const ccomplex_t IT_0135 = IT_0112*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.166666666666667*IT_0013 + -IT_0048 + -IT_0085
       + -IT_0111 + -IT_0136;
    const ccomplex_t IT_0138 = (-0.5)*IT_0013;
    const ccomplex_t IT_0139 = pow(m_sb_1, 2);
    const ccomplex_t IT_0140 = cpow((-2)*s_23 + IT_0006 + IT_0139 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0141 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0021*IT_0023*IT_0025*U_st_00;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = m_b*conj(U_d1)*V_tb*e_em*IT_0025*U_sb_10;
    const ccomplex_t IT_0144 = IT_0024*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0147 = IT_0023*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*(IT_0145 + (-2)*IT_0147);
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = IT_0142*IT_0149;
    const ccomplex_t IT_0151 = IT_0140*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = cpow((-2)*s_23 + IT_0006 + IT_0139 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0154 = m_b*conj(U_d2)*V_tb*e_em*IT_0025*U_sb_10;
    const ccomplex_t IT_0155 = IT_0024*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0158 = IT_0023*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0156 + (-2)*IT_0158);
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0021*IT_0023*IT_0025*U_st_00;
    const ccomplex_t IT_0162 = 0.5*IT_0161;
    const ccomplex_t IT_0163 = IT_0160*IT_0162;
    const ccomplex_t IT_0164 = IT_0153*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = IT_0152 + IT_0165;
    const ccomplex_t IT_0167 = IT_0138 + IT_0166;
    const ccomplex_t IT_0168 = conj(IT_0138) + conj(IT_0166);
    const ccomplex_t IT_0169 = m_b*conj(N_d1)*e_em*IT_0025*U_sb_10;
    const ccomplex_t IT_0170 = IT_0024*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0173 = IT_0016*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0176 = IT_0023*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0171 + 0.333333333333333
      *IT_0174 + -IT_0177);
    const ccomplex_t IT_0179 = 3*IT_0178;
    const ccomplex_t IT_0180 = (-0.333333333333333)*IT_0179;
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = 3*IT_0181;
    const ccomplex_t IT_0183 = (-0.333333333333333)*IT_0182;
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = -IT_0184;
    const ccomplex_t IT_0186 = (-3)*IT_0185;
    const ccomplex_t IT_0187 = (-0.333333333333333)*IT_0186;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = 3*IT_0188;
    const ccomplex_t IT_0190 = 0.166666666666667*IT_0189;
    const ccomplex_t IT_0191 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0192 = IT_0016*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = m_t*N_u1*e_em*IT_0025*U_st_00;
    const ccomplex_t IT_0195 = IT_0036*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + (-0.75)*IT_0196);
    const ccomplex_t IT_0198 = 0.666666666666667*IT_0197;
    const ccomplex_t IT_0199 = IT_0190*IT_0198;
    const ccomplex_t IT_0200 = IT_0014*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = m_t*N_u4*e_em*IT_0025*U_st_00;
    const ccomplex_t IT_0203 = IT_0036*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0206 = IT_0016*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0204 + (
      -1.33333333333333)*IT_0207);
    const ccomplex_t IT_0209 = (-0.75)*IT_0208;
    const ccomplex_t IT_0210 = 0.666666666666667*IT_0209;
    const ccomplex_t IT_0211 = m_b*conj(N_d4)*e_em*IT_0025*U_sb_10;
    const ccomplex_t IT_0212 = IT_0024*IT_0211;
    const ccomplex_t IT_0213 = 1.4142135623731*IT_0212;
    const ccomplex_t IT_0214 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0215 = IT_0023*IT_0214;
    const ccomplex_t IT_0216 = 1.4142135623731*IT_0215;
    const ccomplex_t IT_0217 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0218 = IT_0016*IT_0217;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*(IT_0213 + -IT_0216 +
       0.333333333333333*IT_0219);
    const ccomplex_t IT_0221 = 3*IT_0220;
    const ccomplex_t IT_0222 = 0.333333333333333*IT_0221;
    const ccomplex_t IT_0223 = 3*IT_0222;
    const ccomplex_t IT_0224 = 0.166666666666667*IT_0223;
    const ccomplex_t IT_0225 = IT_0210*IT_0224;
    const ccomplex_t IT_0226 = IT_0112*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0230 = IT_0228*IT_0229;
    const ccomplex_t IT_0231 = IT_0008*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0234 = IT_0016*IT_0233;
    const ccomplex_t IT_0235 = 1.4142135623731*IT_0234;
    const ccomplex_t IT_0236 = m_t*N_u2*e_em*IT_0025*U_st_00;
    const ccomplex_t IT_0237 = IT_0036*IT_0236;
    const ccomplex_t IT_0238 = 1.4142135623731*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*(IT_0235 + (-0.75)*IT_0238);
    const ccomplex_t IT_0240 = 0.666666666666667*IT_0239;
    const ccomplex_t IT_0241 = m_b*conj(N_d2)*e_em*IT_0025*U_sb_10;
    const ccomplex_t IT_0242 = IT_0024*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0245 = IT_0016*IT_0244;
    const ccomplex_t IT_0246 = 1.4142135623731*IT_0245;
    const ccomplex_t IT_0247 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0248 = IT_0023*IT_0247;
    const ccomplex_t IT_0249 = 1.4142135623731*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*(IT_0243 + 0.333333333333333
      *IT_0246 + -IT_0249);
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = 3*IT_0252;
    const ccomplex_t IT_0254 = 0.333333333333333*IT_0253;
    const ccomplex_t IT_0255 = 3*IT_0254;
    const ccomplex_t IT_0256 = 0.333333333333333*IT_0255;
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = -IT_0258;
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = 3*IT_0260;
    const ccomplex_t IT_0262 = 0.333333333333333*IT_0261;
    const ccomplex_t IT_0263 = 3*IT_0262;
    const ccomplex_t IT_0264 = (-0.333333333333333)*IT_0263;
    const ccomplex_t IT_0265 = -IT_0264;
    const ccomplex_t IT_0266 = 3*IT_0265;
    const ccomplex_t IT_0267 = 0.166666666666667*IT_0266;
    const ccomplex_t IT_0268 = IT_0240*IT_0267;
    const ccomplex_t IT_0269 = IT_0049*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*IT_0269;
    const ccomplex_t IT_0271 = m_b*conj(N_d3)*e_em*IT_0025*U_sb_10;
    const ccomplex_t IT_0272 = IT_0024*IT_0271;
    const ccomplex_t IT_0273 = 1.4142135623731*IT_0272;
    const ccomplex_t IT_0274 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0275 = IT_0016*IT_0274;
    const ccomplex_t IT_0276 = 1.4142135623731*IT_0275;
    const ccomplex_t IT_0277 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0278 = IT_0023*IT_0277;
    const ccomplex_t IT_0279 = 1.4142135623731*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*(IT_0273 + 0.333333333333333
      *IT_0276 + -IT_0279);
    const ccomplex_t IT_0281 = -IT_0280;
    const ccomplex_t IT_0282 = -IT_0281;
    const ccomplex_t IT_0283 = -IT_0282;
    const ccomplex_t IT_0284 = -IT_0283;
    const ccomplex_t IT_0285 = 3*IT_0284;
    const ccomplex_t IT_0286 = (-0.333333333333333)*IT_0285;
    const ccomplex_t IT_0287 = -IT_0286;
    const ccomplex_t IT_0288 = 3*IT_0287;
    const ccomplex_t IT_0289 = (-0.333333333333333)*IT_0288;
    const ccomplex_t IT_0290 = -IT_0289;
    const ccomplex_t IT_0291 = -IT_0290;
    const ccomplex_t IT_0292 = -IT_0291;
    const ccomplex_t IT_0293 = -IT_0292;
    const ccomplex_t IT_0294 = -IT_0293;
    const ccomplex_t IT_0295 = 3*IT_0294;
    const ccomplex_t IT_0296 = (-0.333333333333333)*IT_0295;
    const ccomplex_t IT_0297 = (-3)*IT_0296;
    const ccomplex_t IT_0298 = 0.333333333333333*IT_0297;
    const ccomplex_t IT_0299 = 3*IT_0298;
    const ccomplex_t IT_0300 = 0.166666666666667*IT_0299;
    const ccomplex_t IT_0301 = m_t*N_u3*e_em*IT_0025*U_st_00;
    const ccomplex_t IT_0302 = IT_0036*IT_0301;
    const ccomplex_t IT_0303 = 1.4142135623731*IT_0302;
    const ccomplex_t IT_0304 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0305 = IT_0016*IT_0304;
    const ccomplex_t IT_0306 = 1.4142135623731*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*(IT_0303 + (
      -1.33333333333333)*IT_0306);
    const ccomplex_t IT_0308 = (-0.75)*IT_0307;
    const ccomplex_t IT_0309 = (-1.33333333333333)*IT_0308;
    const ccomplex_t IT_0310 = (-0.75)*IT_0309;
    const ccomplex_t IT_0311 = 0.666666666666667*IT_0310;
    const ccomplex_t IT_0312 = IT_0300*IT_0311;
    const ccomplex_t IT_0313 = IT_0086*IT_0312;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*IT_0313;
    const ccomplex_t IT_0315 = -IT_0201 + -IT_0227 + 0.166666666666667*IT_0232
       + -IT_0270 + -IT_0314;
    const ccomplex_t IT_0316 = (-0.5)*IT_0232;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0023*IT_0025*IT_0035*U_sb_00;
    const ccomplex_t IT_0318 = 0.5*IT_0317;
    const ccomplex_t IT_0319 = V_tb*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0320 = IT_0023*IT_0319;
    const ccomplex_t IT_0321 = m_t*V_tb*conj(V_u1)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0322 = IT_0036*IT_0321;
    const ccomplex_t IT_0323 = 1.4142135623731*IT_0322;
    const ccomplex_t IT_0324 = (0 + _Complex_I*1)*(IT_0320 + (-0.5)*IT_0323);
    const ccomplex_t IT_0325 = -IT_0324;
    const ccomplex_t IT_0326 = IT_0318*IT_0325;
    const ccomplex_t IT_0327 = IT_0140*IT_0326;
    const ccomplex_t IT_0328 = (0 + _Complex_I*1)*IT_0327;
    const ccomplex_t IT_0329 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0330 = IT_0023*IT_0329;
    const ccomplex_t IT_0331 = m_t*V_tb*conj(V_u2)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0332 = IT_0036*IT_0331;
    const ccomplex_t IT_0333 = 1.4142135623731*IT_0332;
    const ccomplex_t IT_0334 = (0 + _Complex_I*1)*(IT_0330 + (-0.5)*IT_0333);
    const ccomplex_t IT_0335 = -IT_0334;
    const ccomplex_t IT_0336 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0023*IT_0025*IT_0035*U_sb_00;
    const ccomplex_t IT_0337 = 0.5*IT_0336;
    const ccomplex_t IT_0338 = IT_0335*IT_0337;
    const ccomplex_t IT_0339 = IT_0153*IT_0338;
    const ccomplex_t IT_0340 = (0 + _Complex_I*1)*IT_0339;
    const ccomplex_t IT_0341 = -IT_0328 + -IT_0340;
    const ccomplex_t IT_0342 = IT_0316 + IT_0341;
    const ccomplex_t IT_0343 = conj(IT_0316) + conj(IT_0341);
    const ccomplex_t IT_0344 = (-18)*IT_0001;
    const ccomplex_t IT_0345 = 36*IT_0003;
    const ccomplex_t IT_0346 = IT_0344 + IT_0345;
    const ccomplex_t IT_0347 = m_b*m_t;
    const ccomplex_t IT_0348 = IT_0009*IT_0229;
    const ccomplex_t IT_0349 = IT_0008*IT_0348;
    const ccomplex_t IT_0350 = (0 + _Complex_I*1)*IT_0349;
    const ccomplex_t IT_0351 = m_sG*IT_0350;
    const ccomplex_t IT_0352 = (-0.5)*IT_0351;
    const ccomplex_t IT_0353 = m_b*IT_0166;
    const ccomplex_t IT_0354 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0355 = IT_0149*IT_0325;
    const ccomplex_t IT_0356 = IT_0140*IT_0354*IT_0355;
    const ccomplex_t IT_0357 = m_t*IT_0085;
    const ccomplex_t IT_0358 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0359 = IT_0030*IT_0198;
    const ccomplex_t IT_0360 = IT_0014*IT_0358*IT_0359;
    const ccomplex_t IT_0361 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0362 = IT_0133*IT_0210;
    const ccomplex_t IT_0363 = IT_0112*IT_0361*IT_0362;
    const ccomplex_t IT_0364 = IT_0010*IT_0228;
    const ccomplex_t IT_0365 = IT_0008*IT_0364;
    const ccomplex_t IT_0366 = (0 + _Complex_I*1)*IT_0365;
    const ccomplex_t IT_0367 = m_sG*IT_0366;
    const ccomplex_t IT_0368 = m_t*IT_0013;
    const ccomplex_t IT_0369 = m_t*IT_0048;
    const ccomplex_t IT_0370 = m_t*IT_0111;
    const ccomplex_t IT_0371 = m_t*IT_0136;
    const ccomplex_t IT_0372 = -IT_0357 + IT_0360 + IT_0363 + (
      -0.166666666666667)*IT_0367 + 0.166666666666667*IT_0368 + -IT_0369 + 
      -IT_0370 + -IT_0371;
    const ccomplex_t IT_0373 = 6*IT_0372;
    const ccomplex_t IT_0374 = IT_0142*IT_0318;
    const ccomplex_t IT_0375 = IT_0140*IT_0354*IT_0374;
    const ccomplex_t IT_0376 = (-0.5)*IT_0368 + IT_0375;
    const ccomplex_t IT_0377 = 18*IT_0376;
    const ccomplex_t IT_0378 = 0.5*IT_0367;
    const ccomplex_t IT_0379 = 18*IT_0378;
    const ccomplex_t IT_0380 = m_b*IT_0340;
    const ccomplex_t IT_0381 = -IT_0380;
    const ccomplex_t IT_0382 = (-18)*IT_0381;
    const ccomplex_t IT_0383 = m_b*IT_0328;
    const ccomplex_t IT_0384 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0385 = IT_0162*IT_0337;
    const ccomplex_t IT_0386 = IT_0153*IT_0384*IT_0385;
    const ccomplex_t IT_0387 = -IT_0383 + -IT_0386;
    const ccomplex_t IT_0388 = (-18)*IT_0387;
    const ccomplex_t IT_0389 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0390 = IT_0108*IT_0311;
    const ccomplex_t IT_0391 = IT_0086*IT_0389*IT_0390;
    const ccomplex_t IT_0392 = -IT_0391;
    const ccomplex_t IT_0393 = (-6)*IT_0392;
    const ccomplex_t IT_0394 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0395 = IT_0082*IT_0240;
    const ccomplex_t IT_0396 = IT_0049*IT_0394*IT_0395;
    const ccomplex_t IT_0397 = -IT_0396;
    const ccomplex_t IT_0398 = (-6)*IT_0397;
    const ccomplex_t IT_0399 = IT_0373 + IT_0377 + IT_0379 + IT_0382 + IT_0388
       + IT_0393 + IT_0398;
    const ccomplex_t IT_0400 = 6*conj(IT_0372);
    const ccomplex_t IT_0401 = 18*conj(IT_0376);
    const ccomplex_t IT_0402 = 18*conj(IT_0378);
    const ccomplex_t IT_0403 = (-18)*conj(IT_0381);
    const ccomplex_t IT_0404 = (-18)*conj(IT_0387);
    const ccomplex_t IT_0405 = (-6)*conj(IT_0392);
    const ccomplex_t IT_0406 = (-6)*conj(IT_0397);
    const ccomplex_t IT_0407 = IT_0400 + IT_0401 + IT_0402 + IT_0403 + IT_0404
       + IT_0405 + IT_0406;
    const ccomplex_t IT_0408 = IT_0045*IT_0190;
    const ccomplex_t IT_0409 = IT_0014*IT_0358*IT_0408;
    const ccomplex_t IT_0410 = IT_0123*IT_0224;
    const ccomplex_t IT_0411 = IT_0112*IT_0361*IT_0410;
    const ccomplex_t IT_0412 = m_t*IT_0232;
    const ccomplex_t IT_0413 = m_t*IT_0201;
    const ccomplex_t IT_0414 = m_t*IT_0270;
    const ccomplex_t IT_0415 = m_t*IT_0314;
    const ccomplex_t IT_0416 = m_t*IT_0227;
    const ccomplex_t IT_0417 = IT_0409 + IT_0411 + 0.166666666666667*IT_0412 +
       -IT_0413 + -IT_0414 + -IT_0415 + -IT_0416;
    const ccomplex_t IT_0418 = 6*IT_0381;
    const ccomplex_t IT_0419 = 6*IT_0387;
    const ccomplex_t IT_0420 = 18*IT_0392;
    const ccomplex_t IT_0421 = 18*IT_0397;
    const ccomplex_t IT_0422 = (-18)*IT_0372;
    const ccomplex_t IT_0423 = (-6)*IT_0376;
    const ccomplex_t IT_0424 = (-6)*IT_0378;
    const ccomplex_t IT_0425 = IT_0418 + IT_0419 + IT_0420 + IT_0421 + IT_0422
       + IT_0423 + IT_0424;
    const ccomplex_t IT_0426 = 6*conj(IT_0381);
    const ccomplex_t IT_0427 = 6*conj(IT_0387);
    const ccomplex_t IT_0428 = 18*conj(IT_0392);
    const ccomplex_t IT_0429 = 18*conj(IT_0397);
    const ccomplex_t IT_0430 = (-18)*conj(IT_0372);
    const ccomplex_t IT_0431 = (-6)*conj(IT_0376);
    const ccomplex_t IT_0432 = (-6)*conj(IT_0378);
    const ccomplex_t IT_0433 = IT_0426 + IT_0427 + IT_0428 + IT_0429 + IT_0430
       + IT_0431 + IT_0432;
    const ccomplex_t IT_0434 = IT_0160*IT_0335;
    const ccomplex_t IT_0435 = IT_0153*IT_0384*IT_0434;
    const ccomplex_t IT_0436 = (-0.5)*IT_0412 + -IT_0435;
    const ccomplex_t IT_0437 = 6*IT_0392;
    const ccomplex_t IT_0438 = 6*IT_0397;
    const ccomplex_t IT_0439 = 18*IT_0381;
    const ccomplex_t IT_0440 = 18*IT_0387;
    const ccomplex_t IT_0441 = (-18)*IT_0376;
    const ccomplex_t IT_0442 = (-18)*IT_0378;
    const ccomplex_t IT_0443 = (-6)*IT_0372;
    const ccomplex_t IT_0444 = IT_0437 + IT_0438 + IT_0439 + IT_0440 + IT_0441
       + IT_0442 + IT_0443;
    const ccomplex_t IT_0445 = 6*conj(IT_0392);
    const ccomplex_t IT_0446 = 6*conj(IT_0397);
    const ccomplex_t IT_0447 = 18*conj(IT_0381);
    const ccomplex_t IT_0448 = 18*conj(IT_0387);
    const ccomplex_t IT_0449 = (-18)*conj(IT_0376);
    const ccomplex_t IT_0450 = (-18)*conj(IT_0378);
    const ccomplex_t IT_0451 = (-6)*conj(IT_0372);
    const ccomplex_t IT_0452 = IT_0445 + IT_0446 + IT_0447 + IT_0448 + IT_0449
       + IT_0450 + IT_0451;
    const ccomplex_t IT_0453 = IT_0098*IT_0300;
    const ccomplex_t IT_0454 = IT_0086*IT_0389*IT_0453;
    const ccomplex_t IT_0455 = -IT_0454;
    const ccomplex_t IT_0456 = 0.166666666666667*IT_0351;
    const ccomplex_t IT_0457 = IT_0060*IT_0267;
    const ccomplex_t IT_0458 = IT_0049*IT_0394*IT_0457;
    const ccomplex_t IT_0459 = -IT_0458;
    const ccomplex_t IT_0460 = 6*IT_0376;
    const ccomplex_t IT_0461 = 6*IT_0378;
    const ccomplex_t IT_0462 = 18*IT_0372;
    const ccomplex_t IT_0463 = (-18)*IT_0392;
    const ccomplex_t IT_0464 = (-18)*IT_0397;
    const ccomplex_t IT_0465 = (-6)*IT_0381;
    const ccomplex_t IT_0466 = (-6)*IT_0387;
    const ccomplex_t IT_0467 = IT_0460 + IT_0461 + IT_0462 + IT_0463 + IT_0464
       + IT_0465 + IT_0466;
    const ccomplex_t IT_0468 = 6*conj(IT_0376);
    const ccomplex_t IT_0469 = 6*conj(IT_0378);
    const ccomplex_t IT_0470 = 18*conj(IT_0372);
    const ccomplex_t IT_0471 = (-18)*conj(IT_0392);
    const ccomplex_t IT_0472 = (-18)*conj(IT_0397);
    const ccomplex_t IT_0473 = (-6)*conj(IT_0381);
    const ccomplex_t IT_0474 = (-6)*conj(IT_0387);
    const ccomplex_t IT_0475 = IT_0468 + IT_0469 + IT_0470 + IT_0471 + IT_0472
       + IT_0473 + IT_0474;
    const ccomplex_t IT_0476 = (conj(IT_0352) + conj(IT_0353) + conj(IT_0356))
      *IT_0399 + (IT_0352 + IT_0353 + IT_0356)*IT_0407 + conj(IT_0417)*IT_0425 +
       IT_0417*IT_0433 + conj(IT_0436)*IT_0444 + IT_0436*IT_0452 + (conj(IT_0455
      ) + conj(IT_0456) + conj(IT_0459))*IT_0467 + (IT_0455 + IT_0456 + IT_0459)
      *IT_0475;
    const ccomplex_t IT_0477 = conj(IT_0381) + conj(IT_0387);
    const ccomplex_t IT_0478 = IT_0392 + IT_0397;
    const ccomplex_t IT_0479 = conj(IT_0392) + conj(IT_0397);
    const ccomplex_t IT_0480 = conj(IT_0376) + conj(IT_0378);
    const ccomplex_t IT_0481 = 6*s_34;
    const ccomplex_t IT_0482 = IT_0381 + IT_0387;
    const ccomplex_t IT_0483 = conj(IT_0456) + conj(IT_0459);
    const ccomplex_t IT_0484 = conj(IT_0353) + conj(IT_0356);
    const ccomplex_t IT_0485 = 18*s_34;
    const ccomplex_t IT_0486 = IT_0381*IT_0481;
    const ccomplex_t IT_0487 = IT_0392*IT_0485;
    const ccomplex_t IT_0488 = IT_0376 + IT_0378;
    const ccomplex_t IT_0489 = (-18)*s_34;
    const ccomplex_t IT_0490 = conj(IT_0372)*IT_0489;
    const ccomplex_t IT_0491 = IT_0455*IT_0485;
    const ccomplex_t IT_0492 = IT_0417*IT_0489;
    const ccomplex_t IT_0493 = conj(IT_0417)*IT_0489;
    const ccomplex_t IT_0494 = (-6)*s_34;
    const ccomplex_t IT_0495 = conj(IT_0372)*IT_0494;
    const ccomplex_t IT_0496 = IT_0376*IT_0485;
    const ccomplex_t IT_0497 = IT_0381*IT_0489;
    const ccomplex_t IT_0498 = IT_0392*IT_0494;
    const ccomplex_t IT_0499 = conj(IT_0372)*IT_0481;
    const ccomplex_t IT_0500 = conj(IT_0381)*IT_0489;
    const ccomplex_t IT_0501 = conj(IT_0392)*IT_0494;
    const ccomplex_t IT_0502 = IT_0455*IT_0481;
    const ccomplex_t IT_0503 = IT_0456*IT_0481;
    const ccomplex_t IT_0504 = IT_0502 + IT_0503;
    const ccomplex_t IT_0505 = IT_0352*IT_0485;
    const ccomplex_t IT_0506 = IT_0436*IT_0489;
    const ccomplex_t IT_0507 = IT_0417*IT_0494;
    const ccomplex_t IT_0508 = IT_0459*IT_0481;
    const ccomplex_t IT_0509 = IT_0353*IT_0485;
    const ccomplex_t IT_0510 = conj(IT_0455)*IT_0481;
    const ccomplex_t IT_0511 = conj(IT_0456)*IT_0481;
    const ccomplex_t IT_0512 = IT_0510 + IT_0511;
    const ccomplex_t IT_0513 = conj(IT_0436)*IT_0489;
    const ccomplex_t IT_0514 = conj(IT_0417)*IT_0494;
    const ccomplex_t IT_0515 = IT_0353 + IT_0356;
    const ccomplex_t IT_0516 = conj(IT_0459)*IT_0481;
    const ccomplex_t IT_0517 = conj(IT_0352)*IT_0485;
    const ccomplex_t IT_0518 = IT_0456*IT_0485;
    const ccomplex_t IT_0519 = IT_0436*IT_0494;
    const ccomplex_t IT_0520 = conj(IT_0455)*IT_0485;
    const ccomplex_t IT_0521 = conj(IT_0436)*IT_0494;
    const ccomplex_t IT_0522 = m_t*s_14;
    const ccomplex_t IT_0523 = 6*IT_0315;
    const ccomplex_t IT_0524 = 18*IT_0316;
    const ccomplex_t IT_0525 = 18*IT_0341;
    const ccomplex_t IT_0526 = 6*conj(IT_0315);
    const ccomplex_t IT_0527 = 18*conj(IT_0316);
    const ccomplex_t IT_0528 = 18*conj(IT_0341);
    const ccomplex_t IT_0529 = 6*IT_0455;
    const ccomplex_t IT_0530 = 6*IT_0456;
    const ccomplex_t IT_0531 = 6*IT_0459;
    const ccomplex_t IT_0532 = 18*IT_0352;
    const ccomplex_t IT_0533 = (-18)*IT_0436;
    const ccomplex_t IT_0534 = (-6)*IT_0417;
    const ccomplex_t IT_0535 = 6*conj(IT_0455);
    const ccomplex_t IT_0536 = 6*conj(IT_0456);
    const ccomplex_t IT_0537 = 6*conj(IT_0459);
    const ccomplex_t IT_0538 = 18*conj(IT_0352);
    const ccomplex_t IT_0539 = (-18)*conj(IT_0436);
    const ccomplex_t IT_0540 = (-6)*conj(IT_0417);
    const ccomplex_t IT_0541 = conj(IT_0137)*IT_0425 + IT_0137*IT_0433 +
       IT_0168*IT_0444 + IT_0167*IT_0452 + 6*conj(IT_0315)*(IT_0352 + (-3)
      *IT_0417 + -IT_0436 + 3*IT_0455 + 3*IT_0456 + 3*IT_0459) + 6*IT_0315*(conj
      (IT_0352) + (-3)*conj(IT_0417) + -conj(IT_0436) + 3*conj(IT_0455) + 3*conj
      (IT_0456) + 3*conj(IT_0459)) + IT_0484*(IT_0523 + IT_0524 + IT_0525) +
       IT_0515*(IT_0526 + IT_0527 + IT_0528) + IT_0343*(IT_0529 + IT_0530 +
       IT_0531 + IT_0532 + IT_0533 + IT_0534) + IT_0342*(IT_0535 + IT_0536 +
       IT_0537 + IT_0538 + IT_0539 + IT_0540);
    const ccomplex_t IT_0542 = m_b*s_13;
    const ccomplex_t IT_0543 = (-18)*IT_0138;
    const ccomplex_t IT_0544 = (-18)*IT_0166;
    const ccomplex_t IT_0545 = (-6)*IT_0137;
    const ccomplex_t IT_0546 = IT_0543 + IT_0544 + IT_0545;
    const ccomplex_t IT_0547 = (-18)*conj(IT_0138);
    const ccomplex_t IT_0548 = (-18)*conj(IT_0166);
    const ccomplex_t IT_0549 = (-6)*conj(IT_0137);
    const ccomplex_t IT_0550 = IT_0547 + IT_0548 + IT_0549;
    const ccomplex_t IT_0551 = 6*IT_0417;
    const ccomplex_t IT_0552 = 18*IT_0436;
    const ccomplex_t IT_0553 = (-18)*IT_0352;
    const ccomplex_t IT_0554 = (-6)*IT_0455;
    const ccomplex_t IT_0555 = (-6)*IT_0456;
    const ccomplex_t IT_0556 = (-6)*IT_0459;
    const ccomplex_t IT_0557 = 6*conj(IT_0417);
    const ccomplex_t IT_0558 = 18*conj(IT_0436);
    const ccomplex_t IT_0559 = (-18)*conj(IT_0352);
    const ccomplex_t IT_0560 = (-6)*conj(IT_0455);
    const ccomplex_t IT_0561 = (-6)*conj(IT_0456);
    const ccomplex_t IT_0562 = (-6)*conj(IT_0459);
    const ccomplex_t IT_0563 = IT_0343*IT_0399 + IT_0342*IT_0407 + (-6)*conj
      (IT_0137)*(IT_0352 + (-3)*IT_0417 + -IT_0436 + 3*IT_0455 + 3*IT_0456 + 3
      *IT_0459) + (-6)*IT_0137*(conj(IT_0352) + (-3)*conj(IT_0417) + -conj
      (IT_0436) + 3*conj(IT_0455) + 3*conj(IT_0456) + 3*conj(IT_0459)) + conj
      (IT_0315)*IT_0467 + IT_0315*IT_0475 + IT_0484*IT_0546 + IT_0515*IT_0550 +
       IT_0168*(IT_0551 + IT_0552 + IT_0553 + IT_0554 + IT_0555 + IT_0556) +
       IT_0167*(IT_0557 + IT_0558 + IT_0559 + IT_0560 + IT_0561 + IT_0562);
    const ccomplex_t IT_0564 = m_b*m_t*IT_0000;
    const ccomplex_t IT_0565 = (-18)*IT_0316;
    const ccomplex_t IT_0566 = (-6)*IT_0315;
    const ccomplex_t IT_0567 = (-18)*conj(IT_0316);
    const ccomplex_t IT_0568 = (-6)*conj(IT_0315);
    const ccomplex_t IT_0569 = IT_0564*(conj(IT_0137)*((-18)*IT_0315 + (-6)
      *IT_0316) + IT_0137*((-18)*conj(IT_0315) + (-6)*conj(IT_0316)) + conj
      (IT_0341)*IT_0546 + IT_0341*IT_0550 + IT_0168*(IT_0565 + IT_0566) +
       IT_0167*(IT_0567 + IT_0568));
    const ccomplex_t IT_0570 = IT_0005*(conj(IT_0137)*IT_0167 + IT_0137
      *IT_0168 + conj(IT_0315)*IT_0342 + IT_0315*IT_0343) + (IT_0137*conj
      (IT_0137) + IT_0167*IT_0168 + IT_0315*conj(IT_0315) + IT_0342*IT_0343)
      *IT_0346 + IT_0347*IT_0476 + (conj(IT_0417)*IT_0436 + IT_0417*conj(IT_0436
      ) + conj(IT_0352)*IT_0459 + IT_0352*conj(IT_0459) + IT_0477*IT_0478 +
       IT_0387*IT_0479 + IT_0372*IT_0480)*IT_0481 + (IT_0372*conj(IT_0372) +
       IT_0417*conj(IT_0417) + IT_0436*conj(IT_0436) + IT_0397*IT_0479 + IT_0378
      *IT_0480 + IT_0477*IT_0482 + IT_0459*IT_0483 + IT_0356*IT_0484)*IT_0485 +
       IT_0479*(IT_0486 + IT_0487) + (IT_0372*IT_0479 + IT_0387*IT_0480 + conj
      (IT_0387)*IT_0488)*IT_0489 + IT_0478*IT_0490 + conj(IT_0455)*(IT_0491 +
       IT_0492) + IT_0455*IT_0493 + (conj(IT_0436)*IT_0455 + IT_0436*conj
      (IT_0455) + IT_0372*IT_0477 + IT_0397*IT_0480 + conj(IT_0397)*IT_0488)
      *IT_0494 + IT_0482*IT_0495 + IT_0480*(IT_0496 + IT_0497 + IT_0498) +
       IT_0488*(IT_0499 + IT_0500 + IT_0501) + conj(IT_0352)*(IT_0504 + IT_0505 
      + IT_0506 + IT_0507) + IT_0484*(IT_0504 + IT_0505 + IT_0506 + IT_0507 +
       IT_0508 + IT_0509) + IT_0352*(IT_0512 + IT_0513 + IT_0514) + IT_0515*
      (IT_0512 + IT_0513 + IT_0514 + IT_0516 + IT_0517) + IT_0483*(IT_0491 +
       IT_0492 + IT_0518 + IT_0519) + (IT_0456 + IT_0459)*(IT_0493 + IT_0520 +
       IT_0521) + IT_0522*IT_0541 + IT_0542*IT_0563 + IT_0569;
    return create_ccomplex_return(IT_0570);
}

