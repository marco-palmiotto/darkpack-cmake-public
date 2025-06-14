#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_sb_1_to_t_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_sb_1_to_t_anti_b(
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
    const creal_t m_st_1 = param->m_st_1;
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
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_t, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_st_1, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0010*conj(U_sb_00)*U_st_00;
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
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0022 = cos(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0028 = IT_0010*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = sin(beta);
    const ccomplex_t IT_0031 = cpow(IT_0030, -1);
    const ccomplex_t IT_0032 = IT_0010*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = m_t*conj(N_u2)*e_em*IT_0033*U_st_10;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0026 + 3*IT_0029 + 3
      *IT_0036);
    const ccomplex_t IT_0038 = 0.166666666666667*IT_0037;
    const ccomplex_t IT_0039 = cos(beta);
    const ccomplex_t IT_0040 = cpow(IT_0039, -1);
    const ccomplex_t IT_0041 = IT_0010*IT_0040;
    const ccomplex_t IT_0042 = m_b*N_d2*e_em*IT_0033*conj(U_sb_10);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0046 = IT_0010*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = N_B2*e_em*conj(U_sb_00);
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
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0059 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0060 = IT_0023*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0063 = IT_0010*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = m_t*conj(N_u3)*e_em*IT_0033*U_st_10;
    const ccomplex_t IT_0066 = IT_0032*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0061 + 3*IT_0064 + 3
      *IT_0067);
    const ccomplex_t IT_0069 = 0.166666666666667*IT_0068;
    const ccomplex_t IT_0070 = m_b*N_d3*e_em*IT_0033*conj(U_sb_10);
    const ccomplex_t IT_0071 = IT_0041*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0074 = IT_0023*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0077 = IT_0010*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0072 + 0.333333333333333
      *IT_0075 + -IT_0078);
    const ccomplex_t IT_0080 = 3*IT_0079;
    const ccomplex_t IT_0081 = 0.166666666666667*IT_0080;
    const ccomplex_t IT_0082 = IT_0069*IT_0081;
    const ccomplex_t IT_0083 = IT_0058*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0086 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0087 = IT_0023*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0090 = IT_0010*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = m_t*conj(N_u1)*e_em*IT_0033*U_st_10;
    const ccomplex_t IT_0093 = IT_0032*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0088 + 3*IT_0091 + 3
      *IT_0094);
    const ccomplex_t IT_0096 = 0.166666666666667*IT_0095;
    const ccomplex_t IT_0097 = m_b*N_d1*e_em*IT_0033*conj(U_sb_10);
    const ccomplex_t IT_0098 = IT_0041*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0101 = IT_0010*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0104 = IT_0023*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(IT_0099 + -IT_0102 +
       0.333333333333333*IT_0105);
    const ccomplex_t IT_0107 = 3*IT_0106;
    const ccomplex_t IT_0108 = 0.166666666666667*IT_0107;
    const ccomplex_t IT_0109 = IT_0096*IT_0108;
    const ccomplex_t IT_0110 = IT_0085*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0113 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0114 = IT_0023*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0117 = IT_0010*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = m_t*conj(N_u4)*e_em*IT_0033*U_st_10;
    const ccomplex_t IT_0120 = IT_0032*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0115 + 3*IT_0118 + 3
      *IT_0121);
    const ccomplex_t IT_0123 = 0.166666666666667*IT_0122;
    const ccomplex_t IT_0124 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0125 = IT_0023*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0128 = IT_0010*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = m_b*N_d4*e_em*IT_0033*conj(U_sb_10);
    const ccomplex_t IT_0131 = IT_0041*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0126 + (-3)*IT_0129 + 3
      *IT_0132);
    const ccomplex_t IT_0134 = 0.166666666666667*IT_0133;
    const ccomplex_t IT_0135 = IT_0123*IT_0134;
    const ccomplex_t IT_0136 = IT_0112*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = 0.166666666666667*IT_0008 + -IT_0057 + -IT_0084
       + -IT_0111 + -IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_10);
    const ccomplex_t IT_0141 = IT_0139*IT_0140;
    const ccomplex_t IT_0142 = IT_0003*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = (-0.5)*IT_0143;
    const ccomplex_t IT_0145 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0146 = IT_0023*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = m_t*N_u1*e_em*IT_0033*U_st_00;
    const ccomplex_t IT_0149 = IT_0032*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + (-0.75)*IT_0150);
    const ccomplex_t IT_0152 = 0.666666666666667*IT_0151;
    const ccomplex_t IT_0153 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0154 = IT_0023*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = m_b*conj(N_d1)*e_em*IT_0033*conj(U_sb_00);
    const ccomplex_t IT_0157 = IT_0041*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0155 + 1.5*IT_0158);
    const ccomplex_t IT_0160 = (-0.333333333333333)*IT_0159;
    const ccomplex_t IT_0161 = IT_0152*IT_0160;
    const ccomplex_t IT_0162 = IT_0085*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0165 = IT_0023*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = m_t*N_u2*e_em*IT_0033*U_st_00;
    const ccomplex_t IT_0168 = IT_0032*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*(IT_0166 + (-0.75)*IT_0169);
    const ccomplex_t IT_0171 = 0.666666666666667*IT_0170;
    const ccomplex_t IT_0172 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0173 = IT_0023*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = m_b*conj(N_d2)*e_em*IT_0033*conj(U_sb_00);
    const ccomplex_t IT_0176 = IT_0041*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + 1.5*IT_0177);
    const ccomplex_t IT_0179 = (-0.333333333333333)*IT_0178;
    const ccomplex_t IT_0180 = IT_0171*IT_0179;
    const ccomplex_t IT_0181 = IT_0021*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0184 = IT_0023*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = m_t*N_u3*e_em*IT_0033*U_st_00;
    const ccomplex_t IT_0187 = IT_0032*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(IT_0185 + (-0.75)*IT_0188);
    const ccomplex_t IT_0190 = 0.666666666666667*IT_0189;
    const ccomplex_t IT_0191 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0192 = IT_0023*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = m_b*conj(N_d3)*e_em*IT_0033*conj(U_sb_00);
    const ccomplex_t IT_0195 = IT_0041*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + 1.5*IT_0196);
    const ccomplex_t IT_0198 = (-0.333333333333333)*IT_0197;
    const ccomplex_t IT_0199 = IT_0190*IT_0198;
    const ccomplex_t IT_0200 = IT_0058*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0203 = IT_0023*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = m_t*N_u4*e_em*IT_0033*U_st_00;
    const ccomplex_t IT_0206 = IT_0032*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0204 + (-0.75)*IT_0207);
    const ccomplex_t IT_0209 = 0.666666666666667*IT_0208;
    const ccomplex_t IT_0210 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0211 = IT_0023*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = m_b*conj(N_d4)*e_em*IT_0033*conj(U_sb_00);
    const ccomplex_t IT_0214 = IT_0041*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*(IT_0212 + 1.5*IT_0215);
    const ccomplex_t IT_0217 = (-0.333333333333333)*IT_0216;
    const ccomplex_t IT_0218 = IT_0209*IT_0217;
    const ccomplex_t IT_0219 = IT_0112*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = 0.166666666666667*IT_0143 + -IT_0163 + -IT_0182
       + -IT_0201 + -IT_0220;
    const ccomplex_t IT_0222 = s_34*IT_0002;
    const ccomplex_t IT_0223 = (-6)*IT_0222;
    const ccomplex_t IT_0224 = s_13*s_14;
    const ccomplex_t IT_0225 = 12*IT_0224;
    const ccomplex_t IT_0226 = IT_0223 + IT_0225;
    const ccomplex_t IT_0227 = (-18)*IT_0222;
    const ccomplex_t IT_0228 = 36*IT_0224;
    const ccomplex_t IT_0229 = IT_0227 + IT_0228;
    const ccomplex_t IT_0230 = m_t*s_14;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0010*IT_0033;
    const ccomplex_t IT_0232 = 0.5*IT_0231;
    const ccomplex_t IT_0233 = pow(m_b, 2);
    const ccomplex_t IT_0234 = cpow(IT_0039, 2);
    const ccomplex_t IT_0235 = cpow(IT_0030, 2);
    const ccomplex_t IT_0236 = m_b*mu_h*IT_0030*IT_0033*IT_0040*conj(U_sb_10)
      *U_st_00 + conj(U_sb_00)*((IT_0033*(IT_0000 + -IT_0233) + m_W*(IT_0234 + 
      -IT_0235))*U_st_00 + -m_t*mu_h*IT_0031*IT_0033*IT_0039*U_st_10);
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*V_tb*((-0.707106781186548)
      *e_em*IT_0010*IT_0236 + A_b*IT_0039*conj(U_sb_10)*U_st_00 + -A_t*IT_0030
      *conj(U_sb_00)*U_st_10);
    const ccomplex_t IT_0238 = -IT_0237;
    const ccomplex_t IT_0239 = IT_0232*IT_0238;
    const ccomplex_t IT_0240 = pow(m_W, 2);
    const ccomplex_t IT_0241 = cpow((-2)*s_12 + -IT_0002 + -IT_0015 + IT_0240 
      + -reg_prop, -1);
    const ccomplex_t IT_0242 = IT_0239*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0010*IT_0031*IT_0033*IT_0039;
    const ccomplex_t IT_0245 = 0.5*IT_0244;
    const ccomplex_t IT_0246 = m_t*conj(U_sb_00)*(m_t*IT_0031*IT_0039*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0030*IT_0040*conj(U_sb_00)*U_st_00 + conj
      (U_sb_10)*(mu_h*U_st_00 + m_t*IT_0031*IT_0040*U_st_10));
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0030*conj
      (U_sb_10)*U_st_00 + (-1.4142135623731)*e_em*IT_0010*((-0.5)*IT_0033
      *IT_0246 + m_W*IT_0030*IT_0039*conj(U_sb_00)*U_st_00) + A_t*IT_0039*conj
      (U_sb_00)*U_st_10);
    const ccomplex_t IT_0248 = IT_0245*IT_0247;
    const ccomplex_t IT_0249 = cpow((-2)*s_12 + -IT_0002 + -IT_0015 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0250 = IT_0248*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = IT_0005*IT_0139;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0254 = IT_0003*IT_0252*IT_0253;
    const ccomplex_t IT_0255 = (-0.5)*IT_0011;
    const ccomplex_t IT_0256 = m_t*IT_0255;
    const ccomplex_t IT_0257 = m_t*IT_0012;
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = IT_0256 + IT_0258;
    const ccomplex_t IT_0260 = IT_0018*IT_0259;
    const ccomplex_t IT_0261 = m_t*IT_0008;
    const ccomplex_t IT_0262 = -IT_0243 + -IT_0251 + (-0.5)*IT_0254 + -IT_0260
       + 0.5*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0010*IT_0033;
    const ccomplex_t IT_0264 = (-0.5)*IT_0263;
    const ccomplex_t IT_0265 = IT_0238*IT_0264;
    const ccomplex_t IT_0266 = IT_0241*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0010*IT_0030*IT_0033*IT_0040;
    const ccomplex_t IT_0269 = 0.5*IT_0268;
    const ccomplex_t IT_0270 = IT_0247*IT_0269;
    const ccomplex_t IT_0271 = IT_0249*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = m_t*IT_0143;
    const ccomplex_t IT_0274 = IT_0004*IT_0140;
    const ccomplex_t IT_0275 = IT_0003*IT_0253*IT_0274;
    const ccomplex_t IT_0276 = m_b*IT_0255;
    const ccomplex_t IT_0277 = m_b*IT_0012;
    const ccomplex_t IT_0278 = -IT_0277;
    const ccomplex_t IT_0279 = IT_0276 + IT_0278;
    const ccomplex_t IT_0280 = IT_0018*IT_0279;
    const ccomplex_t IT_0281 = -IT_0267 + -IT_0272 + 0.5*IT_0273 + (-0.5)
      *IT_0275 + IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0283 = IT_0081*IT_0190;
    const ccomplex_t IT_0284 = IT_0058*IT_0282*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0286 = IT_0134*IT_0209;
    const ccomplex_t IT_0287 = IT_0112*IT_0285*IT_0286;
    const ccomplex_t IT_0288 = m_t*IT_0111;
    const ccomplex_t IT_0289 = m_t*IT_0057;
    const ccomplex_t IT_0290 = m_t*IT_0084;
    const ccomplex_t IT_0291 = m_t*IT_0137;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0293 = IT_0108*IT_0152;
    const ccomplex_t IT_0294 = IT_0085*IT_0292*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0296 = IT_0054*IT_0171;
    const ccomplex_t IT_0297 = IT_0021*IT_0295*IT_0296;
    const ccomplex_t IT_0298 = 0.166666666666667*IT_0254 + (-0.166666666666667
      )*IT_0261 + -IT_0284 + -IT_0287 + IT_0288 + IT_0289 + IT_0290 + IT_0291 + 
      -IT_0294 + -IT_0297;
    const ccomplex_t IT_0299 = IT_0020 + 3*IT_0138;
    const ccomplex_t IT_0300 = 3*conj(IT_0138);
    const ccomplex_t IT_0301 = conj(IT_0020) + IT_0300;
    const ccomplex_t IT_0302 = IT_0096*IT_0160;
    const ccomplex_t IT_0303 = IT_0085*IT_0292*IT_0302;
    const ccomplex_t IT_0304 = IT_0123*IT_0217;
    const ccomplex_t IT_0305 = IT_0112*IT_0285*IT_0304;
    const ccomplex_t IT_0306 = m_t*IT_0163;
    const ccomplex_t IT_0307 = m_t*IT_0182;
    const ccomplex_t IT_0308 = m_t*IT_0201;
    const ccomplex_t IT_0309 = m_t*IT_0220;
    const ccomplex_t IT_0310 = IT_0069*IT_0198;
    const ccomplex_t IT_0311 = IT_0058*IT_0282*IT_0310;
    const ccomplex_t IT_0312 = IT_0038*IT_0179;
    const ccomplex_t IT_0313 = IT_0021*IT_0295*IT_0312;
    const ccomplex_t IT_0314 = (-0.166666666666667)*IT_0273 +
       0.166666666666667*IT_0275 + -IT_0303 + -IT_0305 + IT_0306 + IT_0307 +
       IT_0308 + IT_0309 + -IT_0311 + -IT_0313;
    const ccomplex_t IT_0315 = 3*IT_0221;
    const ccomplex_t IT_0316 = IT_0144 + IT_0315;
    const ccomplex_t IT_0317 = 3*conj(IT_0221);
    const ccomplex_t IT_0318 = conj(IT_0144) + IT_0317;
    const ccomplex_t IT_0319 = (18*conj(IT_0020) + 6*conj(IT_0138))*IT_0262 + 
      (18*IT_0020 + 6*IT_0138)*conj(IT_0262) + (18*conj(IT_0144) + 6*conj
      (IT_0221))*IT_0281 + (18*IT_0144 + 6*IT_0221)*conj(IT_0281) + 6*conj
      (IT_0298)*IT_0299 + 6*IT_0298*IT_0301 + 6*conj(IT_0314)*IT_0316 + 6
      *IT_0314*IT_0318;
    const ccomplex_t IT_0320 = 6*s_34;
    const ccomplex_t IT_0321 = 18*s_34;
    const ccomplex_t IT_0322 = m_b*m_t;
    const ccomplex_t IT_0323 = (-18)*conj(IT_0262);
    const ccomplex_t IT_0324 = (-6)*conj(IT_0298);
    const ccomplex_t IT_0325 = IT_0323 + IT_0324;
    const ccomplex_t IT_0326 = (-18)*conj(IT_0298);
    const ccomplex_t IT_0327 = (-6)*conj(IT_0262);
    const ccomplex_t IT_0328 = IT_0326 + IT_0327;
    const ccomplex_t IT_0329 = (-18)*IT_0298;
    const ccomplex_t IT_0330 = (-6)*IT_0298;
    const ccomplex_t IT_0331 = m_b*m_t*IT_0002;
    const ccomplex_t IT_0332 = (-6)*conj(IT_0221);
    const ccomplex_t IT_0333 = (-18)*conj(IT_0144);
    const ccomplex_t IT_0334 = IT_0332 + IT_0333;
    const ccomplex_t IT_0335 = (-6)*IT_0138;
    const ccomplex_t IT_0336 = (-6)*conj(IT_0138);
    const ccomplex_t IT_0337 = m_b*s_13;
    const ccomplex_t IT_0338 = IT_0301*IT_0314 + IT_0299*conj(IT_0314) + (
      -0.166666666666667)*IT_0144*IT_0325 + (-0.166666666666667)*IT_0221*IT_0328
       + (-0.166666666666667)*conj(IT_0221)*IT_0329 + (-0.166666666666667)*conj
      (IT_0144)*IT_0330 + (-0.166666666666667)*IT_0262*IT_0334 + (
      -0.166666666666667)*conj(IT_0281)*((-18)*IT_0020 + IT_0335) + (
      -0.166666666666667)*IT_0281*((-18)*conj(IT_0020) + IT_0336);
    const ccomplex_t IT_0339 = (conj(IT_0020)*IT_0138 + IT_0020*conj(IT_0138) 
      + conj(IT_0144)*IT_0221 + IT_0144*conj(IT_0221))*IT_0226 + (IT_0020*conj
      (IT_0020) + IT_0138*conj(IT_0138) + IT_0144*conj(IT_0144) + IT_0221*conj
      (IT_0221))*IT_0229 + IT_0230*IT_0319 + (conj(IT_0262)*IT_0298 + IT_0262
      *conj(IT_0298) + conj(IT_0281)*IT_0314 + IT_0281*conj(IT_0314))*IT_0320 + 
      (IT_0262*conj(IT_0262) + IT_0281*conj(IT_0281) + IT_0298*conj(IT_0298) +
       IT_0314*conj(IT_0314))*IT_0321 + IT_0322*(IT_0262*((-18)*conj(IT_0281) + 
      (-6)*conj(IT_0314)) + IT_0281*IT_0325 + IT_0314*IT_0328 + conj(IT_0314)
      *IT_0329 + conj(IT_0281)*IT_0330) + IT_0331*(conj(IT_0020)*((-18)*IT_0144 
      + (-6)*IT_0221) + IT_0020*IT_0334 + IT_0318*IT_0335 + IT_0316*IT_0336) + (
      -6)*IT_0337*IT_0338;
    return create_ccomplex_return(IT_0339);
}

