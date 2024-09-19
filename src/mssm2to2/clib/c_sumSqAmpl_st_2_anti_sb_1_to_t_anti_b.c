#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_sb_1_to_t_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_sb_1_to_t_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
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
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_t, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_st_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0010*conj(U_sb_00)*U_st_01;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0010;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = pow(m_sb_1, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0002 + -IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0012*IT_0018;
    const ccomplex_t IT_0020 = (-0.5)*IT_0008 + -IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0022 = cos(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0028 = IT_0010*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = sin(beta);
    const ccomplex_t IT_0031 = cpow(IT_0030, -1);
    const ccomplex_t IT_0032 = IT_0010*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = m_t*conj(N_u3)*e_em*IT_0033*U_st_11;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0026 + 3*IT_0029 + 3
      *IT_0036);
    const ccomplex_t IT_0038 = 0.166666666666667*IT_0037;
    const ccomplex_t IT_0039 = cos(beta);
    const ccomplex_t IT_0040 = cpow(IT_0039, -1);
    const ccomplex_t IT_0041 = IT_0010*IT_0040;
    const ccomplex_t IT_0042 = m_b*N_d3*e_em*IT_0033*conj(U_sb_10);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0046 = IT_0010*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0049 = IT_0023*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0044 + -IT_0047 +
       0.333333333333333*IT_0050);
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = (-3)*IT_0052;
    const ccomplex_t IT_0054 = 0.166666666666667*IT_0053;
    const ccomplex_t IT_0055 = IT_0038*IT_0054;
    const ccomplex_t IT_0056 = IT_0021*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0059 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0060 = IT_0023*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0063 = IT_0010*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = m_t*conj(N_u1)*e_em*IT_0033*U_st_11;
    const ccomplex_t IT_0066 = IT_0032*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0061 + 3*IT_0064 + 3
      *IT_0067);
    const ccomplex_t IT_0069 = 0.166666666666667*IT_0068;
    const ccomplex_t IT_0070 = m_b*N_d1*e_em*IT_0033*conj(U_sb_10);
    const ccomplex_t IT_0071 = IT_0041*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0074 = IT_0023*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0077 = IT_0010*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0072 + 0.333333333333333
      *IT_0075 + -IT_0078);
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = (-3)*IT_0080;
    const ccomplex_t IT_0082 = 0.166666666666667*IT_0081;
    const ccomplex_t IT_0083 = IT_0069*IT_0082;
    const ccomplex_t IT_0084 = IT_0058*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0087 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0088 = IT_0023*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0091 = IT_0010*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = m_t*conj(N_u2)*e_em*IT_0033*U_st_11;
    const ccomplex_t IT_0094 = IT_0032*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0089 + 3*IT_0092 + 3
      *IT_0095);
    const ccomplex_t IT_0097 = 0.166666666666667*IT_0096;
    const ccomplex_t IT_0098 = m_b*N_d2*e_em*IT_0033*conj(U_sb_10);
    const ccomplex_t IT_0099 = IT_0041*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0102 = IT_0010*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0105 = IT_0023*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0100 + -IT_0103 +
       0.333333333333333*IT_0106);
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = (-3)*IT_0108;
    const ccomplex_t IT_0110 = 0.166666666666667*IT_0109;
    const ccomplex_t IT_0111 = IT_0097*IT_0110;
    const ccomplex_t IT_0112 = IT_0086*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0115 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0116 = IT_0023*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0119 = IT_0010*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = m_t*conj(N_u4)*e_em*IT_0033*U_st_11;
    const ccomplex_t IT_0122 = IT_0032*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0117 + 3*IT_0120 + 3
      *IT_0123);
    const ccomplex_t IT_0125 = 0.166666666666667*IT_0124;
    const ccomplex_t IT_0126 = m_b*N_d4*e_em*IT_0033*conj(U_sb_10);
    const ccomplex_t IT_0127 = IT_0041*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0130 = IT_0010*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0133 = IT_0023*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0128 + -IT_0131 +
       0.333333333333333*IT_0134);
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = (-3)*IT_0136;
    const ccomplex_t IT_0138 = 0.166666666666667*IT_0137;
    const ccomplex_t IT_0139 = IT_0125*IT_0138;
    const ccomplex_t IT_0140 = IT_0114*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = 0.166666666666667*IT_0008 + -IT_0057 + -IT_0085
       + -IT_0113 + -IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_10);
    const ccomplex_t IT_0145 = IT_0143*IT_0144;
    const ccomplex_t IT_0146 = IT_0003*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0150 = IT_0023*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = m_b*conj(N_d1)*e_em*IT_0033*conj(U_sb_00);
    const ccomplex_t IT_0153 = IT_0041*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0151 + 1.5*IT_0154);
    const ccomplex_t IT_0156 = (-0.333333333333333)*IT_0155;
    const ccomplex_t IT_0157 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0158 = IT_0023*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = m_t*N_u1*e_em*IT_0033*U_st_01;
    const ccomplex_t IT_0161 = IT_0032*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0159 + (-0.75)*IT_0162);
    const ccomplex_t IT_0164 = 0.666666666666667*IT_0163;
    const ccomplex_t IT_0165 = IT_0156*IT_0164;
    const ccomplex_t IT_0166 = IT_0058*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0169 = IT_0023*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = m_t*N_u2*e_em*IT_0033*U_st_01;
    const ccomplex_t IT_0172 = IT_0032*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0170 + (-0.75)*IT_0173);
    const ccomplex_t IT_0175 = 0.666666666666667*IT_0174;
    const ccomplex_t IT_0176 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0177 = IT_0023*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = m_b*conj(N_d2)*e_em*IT_0033*conj(U_sb_00);
    const ccomplex_t IT_0180 = IT_0041*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*(IT_0178 + 1.5*IT_0181);
    const ccomplex_t IT_0183 = 0.666666666666667*IT_0182;
    const ccomplex_t IT_0184 = 1.5*IT_0183;
    const ccomplex_t IT_0185 = 0.666666666666667*IT_0184;
    const ccomplex_t IT_0186 = 1.5*IT_0185;
    const ccomplex_t IT_0187 = 0.666666666666667*IT_0186;
    const ccomplex_t IT_0188 = 1.5*IT_0187;
    const ccomplex_t IT_0189 = 0.666666666666667*IT_0188;
    const ccomplex_t IT_0190 = 1.5*IT_0189;
    const ccomplex_t IT_0191 = (-0.333333333333333)*IT_0190;
    const ccomplex_t IT_0192 = IT_0175*IT_0191;
    const ccomplex_t IT_0193 = IT_0086*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0196 = IT_0023*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = m_t*N_u3*e_em*IT_0033*U_st_01;
    const ccomplex_t IT_0199 = IT_0032*IT_0198;
    const ccomplex_t IT_0200 = 1.4142135623731*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*(IT_0197 + (-0.75)*IT_0200);
    const ccomplex_t IT_0202 = 0.666666666666667*IT_0201;
    const ccomplex_t IT_0203 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0204 = IT_0023*IT_0203;
    const ccomplex_t IT_0205 = 1.4142135623731*IT_0204;
    const ccomplex_t IT_0206 = m_b*conj(N_d3)*e_em*IT_0033*conj(U_sb_00);
    const ccomplex_t IT_0207 = IT_0041*IT_0206;
    const ccomplex_t IT_0208 = 1.4142135623731*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*(IT_0205 + 1.5*IT_0208);
    const ccomplex_t IT_0210 = 0.666666666666667*IT_0209;
    const ccomplex_t IT_0211 = 1.5*IT_0210;
    const ccomplex_t IT_0212 = (-0.333333333333333)*IT_0211;
    const ccomplex_t IT_0213 = IT_0202*IT_0212;
    const ccomplex_t IT_0214 = IT_0021*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = m_t*N_u4*e_em*IT_0033*U_st_01;
    const ccomplex_t IT_0217 = IT_0032*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0220 = IT_0023*IT_0219;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*(IT_0218 + (
      -1.33333333333333)*IT_0221);
    const ccomplex_t IT_0223 = (-0.75)*IT_0222;
    const ccomplex_t IT_0224 = 0.666666666666667*IT_0223;
    const ccomplex_t IT_0225 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0226 = IT_0023*IT_0225;
    const ccomplex_t IT_0227 = 1.4142135623731*IT_0226;
    const ccomplex_t IT_0228 = m_b*conj(N_d4)*e_em*IT_0033*conj(U_sb_00);
    const ccomplex_t IT_0229 = IT_0041*IT_0228;
    const ccomplex_t IT_0230 = 1.4142135623731*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*(IT_0227 + 1.5*IT_0230);
    const ccomplex_t IT_0232 = 0.666666666666667*IT_0231;
    const ccomplex_t IT_0233 = 1.5*IT_0232;
    const ccomplex_t IT_0234 = (-0.333333333333333)*IT_0233;
    const ccomplex_t IT_0235 = IT_0224*IT_0234;
    const ccomplex_t IT_0236 = IT_0114*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = 0.166666666666667*IT_0147 + -IT_0167 + -IT_0194
       + -IT_0215 + -IT_0237;
    const ccomplex_t IT_0239 = s_34*IT_0002;
    const ccomplex_t IT_0240 = (-6)*IT_0239;
    const ccomplex_t IT_0241 = s_13*s_14;
    const ccomplex_t IT_0242 = 12*IT_0241;
    const ccomplex_t IT_0243 = IT_0240 + IT_0242;
    const ccomplex_t IT_0244 = (-18)*IT_0239;
    const ccomplex_t IT_0245 = 36*IT_0241;
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = m_t*s_14;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0249 = IT_0005*IT_0143;
    const ccomplex_t IT_0250 = IT_0003*IT_0248*IT_0249;
    const ccomplex_t IT_0251 = m_t*IT_0008;
    const ccomplex_t IT_0252 = pow(m_W, 2);
    const ccomplex_t IT_0253 = cpow((-2)*s_12 + -IT_0002 + -IT_0015 + IT_0252 
      + -reg_prop, -1);
    const ccomplex_t IT_0254 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0010*IT_0033;
    const ccomplex_t IT_0255 = 0.5*IT_0254;
    const ccomplex_t IT_0256 = pow(m_b, 2);
    const ccomplex_t IT_0257 = cpow(IT_0039, 2);
    const ccomplex_t IT_0258 = cpow(IT_0030, 2);
    const ccomplex_t IT_0259 = m_b*mu_h*IT_0030*IT_0033*IT_0040*conj(U_sb_10)
      *U_st_01 + conj(U_sb_00)*((IT_0033*(IT_0000 + -IT_0256) + m_W*(IT_0257 + 
      -IT_0258))*U_st_01 + -m_t*mu_h*IT_0031*IT_0033*IT_0039*U_st_11);
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*V_tb*((-0.707106781186548)
      *e_em*IT_0010*IT_0259 + A_b*IT_0039*conj(U_sb_10)*U_st_01 + -A_t*IT_0030
      *conj(U_sb_00)*U_st_11);
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0255*IT_0261;
    const ccomplex_t IT_0263 = IT_0253*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*IT_0263;
    const ccomplex_t IT_0265 = cpow((-2)*s_12 + -IT_0002 + -IT_0015 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0266 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0010*IT_0031*IT_0033*IT_0039;
    const ccomplex_t IT_0267 = 0.5*IT_0266;
    const ccomplex_t IT_0268 = m_t*conj(U_sb_00)*(m_t*IT_0031*IT_0039*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0030*IT_0040*conj(U_sb_00)*U_st_01 + conj
      (U_sb_10)*(mu_h*U_st_01 + m_t*IT_0031*IT_0040*U_st_11));
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0030*conj
      (U_sb_10)*U_st_01 + (-1.4142135623731)*e_em*IT_0010*((-0.5)*IT_0033
      *IT_0268 + m_W*IT_0030*IT_0039*conj(U_sb_00)*U_st_01) + A_t*IT_0039*conj
      (U_sb_00)*U_st_11);
    const ccomplex_t IT_0270 = IT_0267*IT_0269;
    const ccomplex_t IT_0271 = IT_0265*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = (-0.5)*IT_0011;
    const ccomplex_t IT_0274 = m_t*IT_0273;
    const ccomplex_t IT_0275 = m_t*IT_0012;
    const ccomplex_t IT_0276 = -IT_0275;
    const ccomplex_t IT_0277 = IT_0274 + IT_0276;
    const ccomplex_t IT_0278 = IT_0018*IT_0277;
    const ccomplex_t IT_0279 = (-0.5)*IT_0250 + 0.5*IT_0251 + -IT_0264 + 
      -IT_0272 + -IT_0278;
    const ccomplex_t IT_0280 = IT_0004*IT_0144;
    const ccomplex_t IT_0281 = IT_0003*IT_0248*IT_0280;
    const ccomplex_t IT_0282 = m_t*IT_0147;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0010*IT_0033;
    const ccomplex_t IT_0284 = (-0.5)*IT_0283;
    const ccomplex_t IT_0285 = IT_0261*IT_0284;
    const ccomplex_t IT_0286 = IT_0253*IT_0285;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*IT_0286;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0010*IT_0030*IT_0033*IT_0040;
    const ccomplex_t IT_0289 = 0.5*IT_0288;
    const ccomplex_t IT_0290 = IT_0269*IT_0289;
    const ccomplex_t IT_0291 = IT_0265*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = m_b*IT_0273;
    const ccomplex_t IT_0294 = m_b*IT_0012;
    const ccomplex_t IT_0295 = -IT_0294;
    const ccomplex_t IT_0296 = IT_0293 + IT_0295;
    const ccomplex_t IT_0297 = IT_0018*IT_0296;
    const ccomplex_t IT_0298 = (-0.5)*IT_0281 + 0.5*IT_0282 + -IT_0287 + 
      -IT_0292 + IT_0297;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0300 = IT_0082*IT_0164;
    const ccomplex_t IT_0301 = IT_0058*IT_0299*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0303 = IT_0110*IT_0175;
    const ccomplex_t IT_0304 = IT_0086*IT_0302*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0306 = IT_0054*IT_0202;
    const ccomplex_t IT_0307 = IT_0021*IT_0305*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0309 = IT_0138*IT_0224;
    const ccomplex_t IT_0310 = IT_0114*IT_0308*IT_0309;
    const ccomplex_t IT_0311 = m_t*IT_0085;
    const ccomplex_t IT_0312 = m_t*IT_0113;
    const ccomplex_t IT_0313 = m_t*IT_0057;
    const ccomplex_t IT_0314 = m_t*IT_0141;
    const ccomplex_t IT_0315 = 0.166666666666667*IT_0250 + (-0.166666666666667
      )*IT_0251 + -IT_0301 + -IT_0304 + -IT_0307 + -IT_0310 + IT_0311 + IT_0312 
      + IT_0313 + IT_0314;
    const ccomplex_t IT_0316 = IT_0020 + 3*IT_0142;
    const ccomplex_t IT_0317 = 3*conj(IT_0142);
    const ccomplex_t IT_0318 = conj(IT_0020) + IT_0317;
    const ccomplex_t IT_0319 = IT_0038*IT_0212;
    const ccomplex_t IT_0320 = IT_0021*IT_0305*IT_0319;
    const ccomplex_t IT_0321 = IT_0069*IT_0156;
    const ccomplex_t IT_0322 = IT_0058*IT_0299*IT_0321;
    const ccomplex_t IT_0323 = IT_0097*IT_0191;
    const ccomplex_t IT_0324 = IT_0086*IT_0302*IT_0323;
    const ccomplex_t IT_0325 = IT_0125*IT_0234;
    const ccomplex_t IT_0326 = IT_0114*IT_0308*IT_0325;
    const ccomplex_t IT_0327 = m_t*IT_0167;
    const ccomplex_t IT_0328 = m_t*IT_0194;
    const ccomplex_t IT_0329 = m_t*IT_0215;
    const ccomplex_t IT_0330 = m_t*IT_0237;
    const ccomplex_t IT_0331 = 0.166666666666667*IT_0281 + (-0.166666666666667
      )*IT_0282 + -IT_0320 + -IT_0322 + -IT_0324 + -IT_0326 + IT_0327 + IT_0328 
      + IT_0329 + IT_0330;
    const ccomplex_t IT_0332 = 3*IT_0238;
    const ccomplex_t IT_0333 = IT_0148 + IT_0332;
    const ccomplex_t IT_0334 = 3*conj(IT_0238);
    const ccomplex_t IT_0335 = conj(IT_0148) + IT_0334;
    const ccomplex_t IT_0336 = (18*conj(IT_0020) + 6*conj(IT_0142))*IT_0279 + 
      (18*IT_0020 + 6*IT_0142)*conj(IT_0279) + (18*conj(IT_0148) + 6*conj
      (IT_0238))*IT_0298 + (18*IT_0148 + 6*IT_0238)*conj(IT_0298) + 6*conj
      (IT_0315)*IT_0316 + 6*IT_0315*IT_0318 + 6*conj(IT_0331)*IT_0333 + 6
      *IT_0331*IT_0335;
    const ccomplex_t IT_0337 = 6*s_34;
    const ccomplex_t IT_0338 = 18*s_34;
    const ccomplex_t IT_0339 = m_b*m_t;
    const ccomplex_t IT_0340 = (-18)*conj(IT_0279);
    const ccomplex_t IT_0341 = (-6)*conj(IT_0315);
    const ccomplex_t IT_0342 = IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = (-18)*conj(IT_0315);
    const ccomplex_t IT_0344 = (-6)*conj(IT_0279);
    const ccomplex_t IT_0345 = IT_0343 + IT_0344;
    const ccomplex_t IT_0346 = (-18)*IT_0315;
    const ccomplex_t IT_0347 = (-6)*IT_0315;
    const ccomplex_t IT_0348 = m_b*m_t*IT_0002;
    const ccomplex_t IT_0349 = (-6)*conj(IT_0238);
    const ccomplex_t IT_0350 = (-18)*conj(IT_0148);
    const ccomplex_t IT_0351 = IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = (-6)*IT_0142;
    const ccomplex_t IT_0353 = (-6)*conj(IT_0142);
    const ccomplex_t IT_0354 = m_b*s_13;
    const ccomplex_t IT_0355 = IT_0318*IT_0331 + IT_0316*conj(IT_0331) + (
      -0.166666666666667)*IT_0148*IT_0342 + (-0.166666666666667)*IT_0238*IT_0345
       + (-0.166666666666667)*conj(IT_0238)*IT_0346 + (-0.166666666666667)*conj
      (IT_0148)*IT_0347 + (-0.166666666666667)*IT_0279*IT_0351 + (
      -0.166666666666667)*conj(IT_0298)*((-18)*IT_0020 + IT_0352) + (
      -0.166666666666667)*IT_0298*((-18)*conj(IT_0020) + IT_0353);
    const ccomplex_t IT_0356 = (conj(IT_0020)*IT_0142 + IT_0020*conj(IT_0142) 
      + conj(IT_0148)*IT_0238 + IT_0148*conj(IT_0238))*IT_0243 + (IT_0020*conj
      (IT_0020) + IT_0142*conj(IT_0142) + IT_0148*conj(IT_0148) + IT_0238*conj
      (IT_0238))*IT_0246 + IT_0247*IT_0336 + (conj(IT_0279)*IT_0315 + IT_0279
      *conj(IT_0315) + conj(IT_0298)*IT_0331 + IT_0298*conj(IT_0331))*IT_0337 + 
      (IT_0279*conj(IT_0279) + IT_0298*conj(IT_0298) + IT_0315*conj(IT_0315) +
       IT_0331*conj(IT_0331))*IT_0338 + IT_0339*(IT_0279*((-18)*conj(IT_0298) + 
      (-6)*conj(IT_0331)) + IT_0298*IT_0342 + IT_0331*IT_0345 + conj(IT_0331)
      *IT_0346 + conj(IT_0298)*IT_0347) + IT_0348*(conj(IT_0020)*((-18)*IT_0148 
      + (-6)*IT_0238) + IT_0020*IT_0351 + IT_0335*IT_0352 + IT_0333*IT_0353) + (
      -6)*IT_0354*IT_0355;
    return create_ccomplex_return(IT_0356);
}

