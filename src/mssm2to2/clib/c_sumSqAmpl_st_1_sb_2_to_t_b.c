#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_sb_2_to_t_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_sb_2_to_t_b(
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
    const creal_t m_sb_2 = param->m_sb_2;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
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
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = sin(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = sin(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = pow(m_W, -1);
    const ccomplex_t IT_0026 = m_t*N_u1*e_em*IT_0025*U_st_00;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0019 + (-0.75)*IT_0028);
    const ccomplex_t IT_0030 = 0.666666666666667*IT_0029;
    const ccomplex_t IT_0031 = cos(beta);
    const ccomplex_t IT_0032 = cpow(IT_0031, -1);
    const ccomplex_t IT_0033 = IT_0023*IT_0032;
    const ccomplex_t IT_0034 = m_b*conj(N_d1)*e_em*IT_0025*U_sb_11;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0038 = IT_0016*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0041 = IT_0023*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0036 + 0.333333333333333
      *IT_0039 + -IT_0042);
    const ccomplex_t IT_0044 = 3*IT_0043;
    const ccomplex_t IT_0045 = 0.166666666666667*IT_0044;
    const ccomplex_t IT_0046 = IT_0030*IT_0045;
    const ccomplex_t IT_0047 = IT_0014*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0050 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0051 = IT_0016*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = m_t*N_u2*e_em*IT_0025*U_st_00;
    const ccomplex_t IT_0054 = IT_0024*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0052 + (-0.75)*IT_0055);
    const ccomplex_t IT_0057 = 0.666666666666667*IT_0056;
    const ccomplex_t IT_0058 = m_b*conj(N_d2)*e_em*IT_0025*U_sb_11;
    const ccomplex_t IT_0059 = IT_0033*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0062 = IT_0023*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0065 = IT_0016*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0060 + -IT_0063 +
       0.333333333333333*IT_0066);
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = (-3)*IT_0068;
    const ccomplex_t IT_0070 = 0.166666666666667*IT_0069;
    const ccomplex_t IT_0071 = IT_0057*IT_0070;
    const ccomplex_t IT_0072 = IT_0049*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0075 = m_t*N_u3*e_em*IT_0025*U_st_00;
    const ccomplex_t IT_0076 = IT_0024*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0079 = IT_0016*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + (
      -1.33333333333333)*IT_0080);
    const ccomplex_t IT_0082 = (-0.75)*IT_0081;
    const ccomplex_t IT_0083 = (-1.33333333333333)*IT_0082;
    const ccomplex_t IT_0084 = (-0.75)*IT_0083;
    const ccomplex_t IT_0085 = (-1.33333333333333)*IT_0084;
    const ccomplex_t IT_0086 = (-0.75)*IT_0085;
    const ccomplex_t IT_0087 = 0.666666666666667*IT_0086;
    const ccomplex_t IT_0088 = m_b*conj(N_d3)*e_em*IT_0025*U_sb_11;
    const ccomplex_t IT_0089 = IT_0033*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0092 = IT_0016*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0095 = IT_0023*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0090 + 0.333333333333333
      *IT_0093 + -IT_0096);
    const ccomplex_t IT_0098 = 3*IT_0097;
    const ccomplex_t IT_0099 = 0.166666666666667*IT_0098;
    const ccomplex_t IT_0100 = IT_0087*IT_0099;
    const ccomplex_t IT_0101 = IT_0074*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0104 = m_b*conj(N_d4)*e_em*IT_0025*U_sb_11;
    const ccomplex_t IT_0105 = IT_0033*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0108 = IT_0023*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0111 = IT_0016*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0106 + -IT_0109 +
       0.333333333333333*IT_0112);
    const ccomplex_t IT_0114 = 3*IT_0113;
    const ccomplex_t IT_0115 = 0.166666666666667*IT_0114;
    const ccomplex_t IT_0116 = m_t*N_u4*e_em*IT_0025*U_st_00;
    const ccomplex_t IT_0117 = IT_0024*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0120 = IT_0016*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0118 + (
      -1.33333333333333)*IT_0121);
    const ccomplex_t IT_0123 = (-0.75)*IT_0122;
    const ccomplex_t IT_0124 = 0.666666666666667*IT_0123;
    const ccomplex_t IT_0125 = IT_0115*IT_0124;
    const ccomplex_t IT_0126 = IT_0103*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = 0.166666666666667*IT_0013 + -IT_0048 + -IT_0073
       + -IT_0102 + -IT_0127;
    const ccomplex_t IT_0129 = (-0.5)*IT_0013;
    const ccomplex_t IT_0130 = pow(m_sb_2, 2);
    const ccomplex_t IT_0131 = cpow((-2)*s_23 + IT_0006 + IT_0130 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0132 = V_tb*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0133 = IT_0023*IT_0132;
    const ccomplex_t IT_0134 = m_t*V_tb*conj(V_u1)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0135 = IT_0024*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + (-0.5)*IT_0136);
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0021*IT_0023*IT_0025*U_sb_01;
    const ccomplex_t IT_0140 = 0.5*IT_0139;
    const ccomplex_t IT_0141 = IT_0138*IT_0140;
    const ccomplex_t IT_0142 = IT_0131*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = cpow((-2)*s_23 + IT_0006 + IT_0130 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0145 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0146 = IT_0023*IT_0145;
    const ccomplex_t IT_0147 = m_t*V_tb*conj(V_u2)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0148 = IT_0024*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0146 + (-0.5)*IT_0149);
    const ccomplex_t IT_0151 = -IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0021*IT_0023*IT_0025*U_sb_01;
    const ccomplex_t IT_0153 = 0.5*IT_0152;
    const ccomplex_t IT_0154 = IT_0151*IT_0153;
    const ccomplex_t IT_0155 = IT_0144*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = -IT_0143 + -IT_0156;
    const ccomplex_t IT_0158 = IT_0129 + IT_0157;
    const ccomplex_t IT_0159 = conj(IT_0129) + conj(IT_0157);
    const ccomplex_t IT_0160 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0162 = IT_0160*IT_0161;
    const ccomplex_t IT_0163 = IT_0008*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0166 = IT_0023*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = m_t*conj(N_u1)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0169 = IT_0024*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0172 = IT_0016*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0167 + IT_0170 +
       0.333333333333333*IT_0173);
    const ccomplex_t IT_0175 = 3*IT_0174;
    const ccomplex_t IT_0176 = 0.166666666666667*IT_0175;
    const ccomplex_t IT_0177 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0178 = IT_0016*IT_0177;
    const ccomplex_t IT_0179 = 1.4142135623731*IT_0178;
    const ccomplex_t IT_0180 = m_b*N_d1*e_em*IT_0025*U_sb_01;
    const ccomplex_t IT_0181 = IT_0033*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0179 + 1.5*IT_0182);
    const ccomplex_t IT_0184 = (-0.333333333333333)*IT_0183;
    const ccomplex_t IT_0185 = IT_0176*IT_0184;
    const ccomplex_t IT_0186 = IT_0014*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0189 = IT_0016*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0192 = IT_0023*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = m_t*conj(N_u2)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0195 = IT_0024*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0190 + 3*IT_0193 + 3
      *IT_0196);
    const ccomplex_t IT_0198 = 0.166666666666667*IT_0197;
    const ccomplex_t IT_0199 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0200 = IT_0016*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = m_b*N_d2*e_em*IT_0025*U_sb_01;
    const ccomplex_t IT_0203 = IT_0033*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*(IT_0201 + 1.5*IT_0204);
    const ccomplex_t IT_0206 = 0.666666666666667*IT_0205;
    const ccomplex_t IT_0207 = 1.5*IT_0206;
    const ccomplex_t IT_0208 = (-0.333333333333333)*IT_0207;
    const ccomplex_t IT_0209 = IT_0198*IT_0208;
    const ccomplex_t IT_0210 = IT_0049*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0213 = IT_0016*IT_0212;
    const ccomplex_t IT_0214 = 1.4142135623731*IT_0213;
    const ccomplex_t IT_0215 = m_b*N_d3*e_em*IT_0025*U_sb_01;
    const ccomplex_t IT_0216 = IT_0033*IT_0215;
    const ccomplex_t IT_0217 = 1.4142135623731*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*(IT_0214 + 1.5*IT_0217);
    const ccomplex_t IT_0219 = 0.666666666666667*IT_0218;
    const ccomplex_t IT_0220 = 1.5*IT_0219;
    const ccomplex_t IT_0221 = 0.666666666666667*IT_0220;
    const ccomplex_t IT_0222 = 1.5*IT_0221;
    const ccomplex_t IT_0223 = 0.666666666666667*IT_0222;
    const ccomplex_t IT_0224 = 1.5*IT_0223;
    const ccomplex_t IT_0225 = (-0.333333333333333)*IT_0224;
    const ccomplex_t IT_0226 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0227 = IT_0023*IT_0226;
    const ccomplex_t IT_0228 = 1.4142135623731*IT_0227;
    const ccomplex_t IT_0229 = m_t*conj(N_u3)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0230 = IT_0024*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0233 = IT_0016*IT_0232;
    const ccomplex_t IT_0234 = 1.4142135623731*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*(IT_0228 + IT_0231 +
       0.333333333333333*IT_0234);
    const ccomplex_t IT_0236 = 3*IT_0235;
    const ccomplex_t IT_0237 = 0.166666666666667*IT_0236;
    const ccomplex_t IT_0238 = IT_0225*IT_0237;
    const ccomplex_t IT_0239 = IT_0074*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*IT_0239;
    const ccomplex_t IT_0241 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0242 = IT_0023*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = m_t*conj(N_u4)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0245 = IT_0024*IT_0244;
    const ccomplex_t IT_0246 = 1.4142135623731*IT_0245;
    const ccomplex_t IT_0247 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0248 = IT_0016*IT_0247;
    const ccomplex_t IT_0249 = 1.4142135623731*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*(IT_0243 + IT_0246 +
       0.333333333333333*IT_0249);
    const ccomplex_t IT_0251 = 3*IT_0250;
    const ccomplex_t IT_0252 = 0.166666666666667*IT_0251;
    const ccomplex_t IT_0253 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0254 = IT_0016*IT_0253;
    const ccomplex_t IT_0255 = 1.4142135623731*IT_0254;
    const ccomplex_t IT_0256 = m_b*N_d4*e_em*IT_0025*U_sb_01;
    const ccomplex_t IT_0257 = IT_0033*IT_0256;
    const ccomplex_t IT_0258 = 1.4142135623731*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*(IT_0255 + 1.5*IT_0258);
    const ccomplex_t IT_0260 = 0.666666666666667*IT_0259;
    const ccomplex_t IT_0261 = 1.5*IT_0260;
    const ccomplex_t IT_0262 = (-0.333333333333333)*IT_0261;
    const ccomplex_t IT_0263 = IT_0252*IT_0262;
    const ccomplex_t IT_0264 = IT_0103*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = 0.166666666666667*IT_0164 + -IT_0187 + -IT_0211
       + -IT_0240 + -IT_0265;
    const ccomplex_t IT_0267 = (-0.5)*IT_0164;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0023*IT_0025*IT_0032*U_st_00;
    const ccomplex_t IT_0269 = 0.5*IT_0268;
    const ccomplex_t IT_0270 = m_b*conj(U_d1)*V_tb*e_em*IT_0025*U_sb_11;
    const ccomplex_t IT_0271 = IT_0033*IT_0270;
    const ccomplex_t IT_0272 = 1.4142135623731*IT_0271;
    const ccomplex_t IT_0273 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0274 = IT_0023*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*(IT_0272 + (-2)*IT_0274);
    const ccomplex_t IT_0276 = (-0.5)*IT_0275;
    const ccomplex_t IT_0277 = IT_0269*IT_0276;
    const ccomplex_t IT_0278 = IT_0131*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0023*IT_0025*IT_0032*U_st_00;
    const ccomplex_t IT_0281 = 0.5*IT_0280;
    const ccomplex_t IT_0282 = m_b*conj(U_d2)*V_tb*e_em*IT_0025*U_sb_11;
    const ccomplex_t IT_0283 = IT_0033*IT_0282;
    const ccomplex_t IT_0284 = 1.4142135623731*IT_0283;
    const ccomplex_t IT_0285 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0286 = IT_0023*IT_0285;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*(IT_0284 + (-2)*IT_0286);
    const ccomplex_t IT_0288 = (-0.5)*IT_0287;
    const ccomplex_t IT_0289 = IT_0281*IT_0288;
    const ccomplex_t IT_0290 = IT_0144*IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*IT_0290;
    const ccomplex_t IT_0292 = IT_0279 + IT_0291;
    const ccomplex_t IT_0293 = IT_0267 + IT_0292;
    const ccomplex_t IT_0294 = conj(IT_0267) + conj(IT_0292);
    const ccomplex_t IT_0295 = (-18)*IT_0001;
    const ccomplex_t IT_0296 = 36*IT_0003;
    const ccomplex_t IT_0297 = IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = m_b*m_t;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0300 = IT_0099*IT_0237;
    const ccomplex_t IT_0301 = IT_0074*IT_0299*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0303 = IT_0045*IT_0176;
    const ccomplex_t IT_0304 = IT_0014*IT_0302*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0306 = IT_0115*IT_0252;
    const ccomplex_t IT_0307 = IT_0103*IT_0305*IT_0306;
    const ccomplex_t IT_0308 = m_t*IT_0013;
    const ccomplex_t IT_0309 = m_t*IT_0048;
    const ccomplex_t IT_0310 = m_t*IT_0073;
    const ccomplex_t IT_0311 = m_t*IT_0102;
    const ccomplex_t IT_0312 = m_t*IT_0127;
    const ccomplex_t IT_0313 = IT_0301 + IT_0304 + IT_0307 + 0.166666666666667
      *IT_0308 + -IT_0309 + -IT_0310 + -IT_0311 + -IT_0312;
    const ccomplex_t IT_0314 = IT_0143 + IT_0156;
    const ccomplex_t IT_0315 = m_b*IT_0314;
    const ccomplex_t IT_0316 = -IT_0315;
    const ccomplex_t IT_0317 = 6*IT_0316;
    const ccomplex_t IT_0318 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0319 = IT_0057*IT_0208;
    const ccomplex_t IT_0320 = IT_0049*IT_0318*IT_0319;
    const ccomplex_t IT_0321 = -IT_0320;
    const ccomplex_t IT_0322 = 18*IT_0321;
    const ccomplex_t IT_0323 = m_t*IT_0164;
    const ccomplex_t IT_0324 = IT_0030*IT_0184;
    const ccomplex_t IT_0325 = IT_0014*IT_0302*IT_0324;
    const ccomplex_t IT_0326 = IT_0087*IT_0225;
    const ccomplex_t IT_0327 = IT_0074*IT_0299*IT_0326;
    const ccomplex_t IT_0328 = IT_0124*IT_0262;
    const ccomplex_t IT_0329 = IT_0103*IT_0305*IT_0328;
    const ccomplex_t IT_0330 = IT_0009*IT_0161;
    const ccomplex_t IT_0331 = IT_0008*IT_0330;
    const ccomplex_t IT_0332 = (0 + _Complex_I*1)*IT_0331;
    const ccomplex_t IT_0333 = m_sG*IT_0332;
    const ccomplex_t IT_0334 = m_t*IT_0187;
    const ccomplex_t IT_0335 = m_t*IT_0211;
    const ccomplex_t IT_0336 = m_t*IT_0240;
    const ccomplex_t IT_0337 = m_t*IT_0265;
    const ccomplex_t IT_0338 = 0.166666666666667*IT_0323 + IT_0325 + IT_0327 +
       IT_0329 + (-0.166666666666667)*IT_0333 + -IT_0334 + -IT_0335 + -IT_0336 +
       -IT_0337;
    const ccomplex_t IT_0339 = (-18)*IT_0338;
    const ccomplex_t IT_0340 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0341 = IT_0140*IT_0269;
    const ccomplex_t IT_0342 = IT_0131*IT_0340*IT_0341;
    const ccomplex_t IT_0343 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0344 = IT_0153*IT_0281;
    const ccomplex_t IT_0345 = IT_0144*IT_0343*IT_0344;
    const ccomplex_t IT_0346 = (-0.5)*IT_0323 + IT_0342 + IT_0345;
    const ccomplex_t IT_0347 = (-6)*IT_0346;
    const ccomplex_t IT_0348 = 0.5*IT_0333;
    const ccomplex_t IT_0349 = (-6)*IT_0348;
    const ccomplex_t IT_0350 = IT_0317 + IT_0322 + IT_0339 + IT_0347 + IT_0349;
    const ccomplex_t IT_0351 = 6*conj(IT_0316);
    const ccomplex_t IT_0352 = 18*conj(IT_0321);
    const ccomplex_t IT_0353 = (-18)*conj(IT_0338);
    const ccomplex_t IT_0354 = (-6)*conj(IT_0346);
    const ccomplex_t IT_0355 = (-6)*conj(IT_0348);
    const ccomplex_t IT_0356 = IT_0351 + IT_0352 + IT_0353 + IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = IT_0010*IT_0160;
    const ccomplex_t IT_0358 = IT_0008*IT_0357;
    const ccomplex_t IT_0359 = (0 + _Complex_I*1)*IT_0358;
    const ccomplex_t IT_0360 = m_sG*IT_0359;
    const ccomplex_t IT_0361 = m_b*IT_0279;
    const ccomplex_t IT_0362 = m_b*IT_0291;
    const ccomplex_t IT_0363 = (-0.5)*IT_0360 + IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = 6*IT_0338;
    const ccomplex_t IT_0365 = 18*IT_0346;
    const ccomplex_t IT_0366 = 18*IT_0348;
    const ccomplex_t IT_0367 = (-18)*IT_0316;
    const ccomplex_t IT_0368 = (-6)*IT_0321;
    const ccomplex_t IT_0369 = IT_0364 + IT_0365 + IT_0366 + IT_0367 + IT_0368;
    const ccomplex_t IT_0370 = 6*conj(IT_0338);
    const ccomplex_t IT_0371 = 18*conj(IT_0346);
    const ccomplex_t IT_0372 = 18*conj(IT_0348);
    const ccomplex_t IT_0373 = (-18)*conj(IT_0316);
    const ccomplex_t IT_0374 = (-6)*conj(IT_0321);
    const ccomplex_t IT_0375 = IT_0370 + IT_0371 + IT_0372 + IT_0373 + IT_0374;
    const ccomplex_t IT_0376 = IT_0070*IT_0198;
    const ccomplex_t IT_0377 = IT_0049*IT_0318*IT_0376;
    const ccomplex_t IT_0378 = 0.166666666666667*IT_0360 + -IT_0377;
    const ccomplex_t IT_0379 = 6*IT_0346;
    const ccomplex_t IT_0380 = 6*IT_0348;
    const ccomplex_t IT_0381 = 18*IT_0338;
    const ccomplex_t IT_0382 = (-18)*IT_0321;
    const ccomplex_t IT_0383 = (-6)*IT_0316;
    const ccomplex_t IT_0384 = IT_0379 + IT_0380 + IT_0381 + IT_0382 + IT_0383;
    const ccomplex_t IT_0385 = 6*conj(IT_0346);
    const ccomplex_t IT_0386 = 6*conj(IT_0348);
    const ccomplex_t IT_0387 = 18*conj(IT_0338);
    const ccomplex_t IT_0388 = (-18)*conj(IT_0321);
    const ccomplex_t IT_0389 = (-6)*conj(IT_0316);
    const ccomplex_t IT_0390 = IT_0385 + IT_0386 + IT_0387 + IT_0388 + IT_0389;
    const ccomplex_t IT_0391 = IT_0138*IT_0276;
    const ccomplex_t IT_0392 = IT_0131*IT_0340*IT_0391;
    const ccomplex_t IT_0393 = IT_0151*IT_0288;
    const ccomplex_t IT_0394 = IT_0144*IT_0343*IT_0393;
    const ccomplex_t IT_0395 = -IT_0392 + -IT_0394;
    const ccomplex_t IT_0396 = (-0.5)*IT_0308;
    const ccomplex_t IT_0397 = conj(IT_0395) + conj(IT_0396);
    const ccomplex_t IT_0398 = 6*IT_0321;
    const ccomplex_t IT_0399 = 18*IT_0316;
    const ccomplex_t IT_0400 = (-18)*IT_0346;
    const ccomplex_t IT_0401 = (-18)*IT_0348;
    const ccomplex_t IT_0402 = (-6)*IT_0338;
    const ccomplex_t IT_0403 = IT_0398 + IT_0399 + IT_0400 + IT_0401 + IT_0402;
    const ccomplex_t IT_0404 = IT_0395 + IT_0396;
    const ccomplex_t IT_0405 = 6*conj(IT_0321);
    const ccomplex_t IT_0406 = 18*conj(IT_0316);
    const ccomplex_t IT_0407 = (-18)*conj(IT_0346);
    const ccomplex_t IT_0408 = (-18)*conj(IT_0348);
    const ccomplex_t IT_0409 = (-6)*conj(IT_0338);
    const ccomplex_t IT_0410 = IT_0405 + IT_0406 + IT_0407 + IT_0408 + IT_0409;
    const ccomplex_t IT_0411 = IT_0346 + IT_0348;
    const ccomplex_t IT_0412 = 6*s_34;
    const ccomplex_t IT_0413 = conj(IT_0346) + conj(IT_0348);
    const ccomplex_t IT_0414 = 18*s_34;
    const ccomplex_t IT_0415 = IT_0338*IT_0412;
    const ccomplex_t IT_0416 = IT_0346*IT_0414;
    const ccomplex_t IT_0417 = (-18)*s_34;
    const ccomplex_t IT_0418 = conj(IT_0316)*IT_0417;
    const ccomplex_t IT_0419 = IT_0313*IT_0412;
    const ccomplex_t IT_0420 = IT_0395*IT_0414;
    const ccomplex_t IT_0421 = IT_0363*IT_0417;
    const ccomplex_t IT_0422 = conj(IT_0313)*IT_0412;
    const ccomplex_t IT_0423 = conj(IT_0363)*IT_0417;
    const ccomplex_t IT_0424 = (-6)*s_34;
    const ccomplex_t IT_0425 = m_t*s_14;
    const ccomplex_t IT_0426 = -IT_0396;
    const ccomplex_t IT_0427 = (-3)*IT_0313 + IT_0363 + 3*IT_0378 + -IT_0395 +
       IT_0426;
    const ccomplex_t IT_0428 = -conj(IT_0396);
    const ccomplex_t IT_0429 = -conj(IT_0395);
    const ccomplex_t IT_0430 = (-3)*conj(IT_0313);
    const ccomplex_t IT_0431 = conj(IT_0363) + 3*conj(IT_0378) + IT_0428 +
       IT_0429 + IT_0430;
    const ccomplex_t IT_0432 = 6*IT_0378;
    const ccomplex_t IT_0433 = 18*IT_0363;
    const ccomplex_t IT_0434 = (-18)*IT_0395;
    const ccomplex_t IT_0435 = (-18)*IT_0396;
    const ccomplex_t IT_0436 = (-6)*IT_0313;
    const ccomplex_t IT_0437 = 6*conj(IT_0378);
    const ccomplex_t IT_0438 = 18*conj(IT_0363);
    const ccomplex_t IT_0439 = (-18)*conj(IT_0395);
    const ccomplex_t IT_0440 = (-18)*conj(IT_0396);
    const ccomplex_t IT_0441 = (-6)*conj(IT_0313);
    const ccomplex_t IT_0442 = 0.166666666666667*conj(IT_0266)*IT_0350 +
       0.166666666666667*IT_0266*IT_0356 + 0.166666666666667*IT_0294*IT_0403 +
       0.166666666666667*IT_0293*IT_0410 + conj(IT_0128)*IT_0427 + IT_0128
      *IT_0431 + 0.166666666666667*IT_0159*(IT_0432 + IT_0433 + IT_0434 +
       IT_0435 + IT_0436) + 0.166666666666667*IT_0158*(IT_0437 + IT_0438 +
       IT_0439 + IT_0440 + IT_0441);
    const ccomplex_t IT_0443 = m_b*s_13;
    const ccomplex_t IT_0444 = 6*IT_0313;
    const ccomplex_t IT_0445 = 18*IT_0395;
    const ccomplex_t IT_0446 = 18*IT_0396;
    const ccomplex_t IT_0447 = (-18)*IT_0363;
    const ccomplex_t IT_0448 = (-6)*IT_0378;
    const ccomplex_t IT_0449 = 6*conj(IT_0313);
    const ccomplex_t IT_0450 = 18*conj(IT_0395);
    const ccomplex_t IT_0451 = 18*conj(IT_0396);
    const ccomplex_t IT_0452 = (-18)*conj(IT_0363);
    const ccomplex_t IT_0453 = (-6)*conj(IT_0378);
    const ccomplex_t IT_0454 = (-0.166666666666667)*IT_0159*IT_0369 + (
      -0.166666666666667)*IT_0158*IT_0375 + (-0.166666666666667)*conj(IT_0128)
      *IT_0384 + (-0.166666666666667)*IT_0128*IT_0390 + conj(IT_0266)*IT_0427 +
       IT_0266*IT_0431 + (-0.166666666666667)*IT_0294*(IT_0444 + IT_0445 +
       IT_0446 + IT_0447 + IT_0448) + (-0.166666666666667)*IT_0293*(IT_0449 +
       IT_0450 + IT_0451 + IT_0452 + IT_0453);
    const ccomplex_t IT_0455 = m_b*m_t*IT_0000;
    const ccomplex_t IT_0456 = (-18)*IT_0129;
    const ccomplex_t IT_0457 = (-6)*IT_0128;
    const ccomplex_t IT_0458 = (-18)*conj(IT_0129);
    const ccomplex_t IT_0459 = (-6)*conj(IT_0128);
    const ccomplex_t IT_0460 = ((-18)*conj(IT_0128) + (-6)*conj(IT_0129) + (-6
      )*conj(IT_0157))*IT_0266 + ((-18)*IT_0128 + (-6)*IT_0129 + (-6)*IT_0157)
      *conj(IT_0266) + IT_0294*((-18)*IT_0157 + IT_0456 + IT_0457) + IT_0293*((
      -18)*conj(IT_0157) + IT_0458 + IT_0459);
    const ccomplex_t IT_0461 = IT_0005*(conj(IT_0128)*IT_0158 + IT_0128
      *IT_0159 + conj(IT_0266)*IT_0293 + IT_0266*IT_0294) + (IT_0128*conj
      (IT_0128) + IT_0158*IT_0159 + IT_0266*conj(IT_0266) + IT_0293*IT_0294)
      *IT_0297 + IT_0298*(conj(IT_0313)*IT_0350 + IT_0313*IT_0356 + conj(IT_0363
      )*IT_0369 + IT_0363*IT_0375 + conj(IT_0378)*IT_0384 + IT_0378*IT_0390 +
       IT_0397*IT_0403 + IT_0404*IT_0410) + (conj(IT_0316)*IT_0321 + IT_0316
      *conj(IT_0321) + conj(IT_0363)*IT_0378 + IT_0363*conj(IT_0378) + conj
      (IT_0338)*IT_0411)*IT_0412 + (IT_0313*conj(IT_0313) + IT_0316*conj(IT_0316
      ) + IT_0321*conj(IT_0321) + IT_0338*conj(IT_0338) + IT_0363*conj(IT_0363) 
      + IT_0378*conj(IT_0378) + IT_0396*IT_0397 + IT_0348*IT_0413)*IT_0414 +
       IT_0413*(IT_0415 + IT_0416) + (conj(IT_0321)*IT_0338 + IT_0321*conj
      (IT_0338) + conj(IT_0313)*IT_0378 + IT_0313*conj(IT_0378) + IT_0316
      *IT_0413)*IT_0417 + IT_0411*IT_0418 + IT_0397*(IT_0419 + IT_0420 + IT_0421
      ) + IT_0404*(IT_0422 + IT_0423) + (conj(IT_0316)*IT_0338 + IT_0316*conj
      (IT_0338) + conj(IT_0313)*IT_0363 + IT_0313*conj(IT_0363) + IT_0378
      *IT_0397 + conj(IT_0378)*IT_0404 + conj(IT_0321)*IT_0411 + IT_0321*IT_0413
      )*IT_0424 + 6*IT_0425*IT_0442 + (-6)*IT_0443*IT_0454 + IT_0455*IT_0460;
    return create_ccomplex_return(IT_0461);
}

