#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_sb_2_to_t_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_sb_2_to_t_b(
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
    const creal_t m_st_2 = param->m_st_2;
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
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_st_2, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-6)*IT_0001;
    const ccomplex_t IT_0003 = s_13*s_14;
    const ccomplex_t IT_0004 = 12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_t, 2);
    const ccomplex_t IT_0007 = pow(m_sG, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -IT_0007 +
       reg_prop, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = m_b*conj(N_d4)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0025 = IT_0018*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = cos(theta_W);
    const ccomplex_t IT_0028 = cpow(IT_0027, -1);
    const ccomplex_t IT_0029 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0023 + -IT_0026 +
       0.333333333333333*IT_0031);
    const ccomplex_t IT_0033 = 3*IT_0032;
    const ccomplex_t IT_0034 = 0.166666666666667*IT_0033;
    const ccomplex_t IT_0035 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0036 = IT_0028*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = sin(beta);
    const ccomplex_t IT_0039 = cpow(IT_0038, -1);
    const ccomplex_t IT_0040 = IT_0018*IT_0039;
    const ccomplex_t IT_0041 = m_t*N_u4*e_em*IT_0020*U_st_01;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0037 + (-0.75)*IT_0043);
    const ccomplex_t IT_0045 = (-1.33333333333333)*IT_0044;
    const ccomplex_t IT_0046 = (-0.75)*IT_0045;
    const ccomplex_t IT_0047 = 0.666666666666667*IT_0046;
    const ccomplex_t IT_0048 = IT_0034*IT_0047;
    const ccomplex_t IT_0049 = IT_0014*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0052 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0053 = IT_0028*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = m_t*N_u1*e_em*IT_0020*U_st_01;
    const ccomplex_t IT_0056 = IT_0040*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0054 + (-0.75)*IT_0057);
    const ccomplex_t IT_0059 = 0.666666666666667*IT_0058;
    const ccomplex_t IT_0060 = m_b*conj(N_d1)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0061 = IT_0019*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0064 = IT_0028*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0067 = IT_0018*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0062 + 0.333333333333333
      *IT_0065 + -IT_0068);
    const ccomplex_t IT_0070 = 3*IT_0069;
    const ccomplex_t IT_0071 = 0.166666666666667*IT_0070;
    const ccomplex_t IT_0072 = IT_0059*IT_0071;
    const ccomplex_t IT_0073 = IT_0051*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0076 = m_b*conj(N_d2)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0077 = IT_0019*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0080 = IT_0018*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0083 = IT_0028*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0078 + -IT_0081 +
       0.333333333333333*IT_0084);
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = (-3)*IT_0086;
    const ccomplex_t IT_0088 = 0.166666666666667*IT_0087;
    const ccomplex_t IT_0089 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0090 = IT_0028*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = m_t*N_u2*e_em*IT_0020*U_st_01;
    const ccomplex_t IT_0093 = IT_0040*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0091 + (-0.75)*IT_0094);
    const ccomplex_t IT_0096 = 0.666666666666667*IT_0095;
    const ccomplex_t IT_0097 = IT_0088*IT_0096;
    const ccomplex_t IT_0098 = IT_0075*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0101 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0102 = IT_0028*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = m_t*N_u3*e_em*IT_0020*U_st_01;
    const ccomplex_t IT_0105 = IT_0040*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0103 + (-0.75)*IT_0106);
    const ccomplex_t IT_0108 = 0.666666666666667*IT_0107;
    const ccomplex_t IT_0109 = m_b*conj(N_d3)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0110 = IT_0019*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0113 = IT_0028*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0116 = IT_0018*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0111 + 0.333333333333333
      *IT_0114 + -IT_0117);
    const ccomplex_t IT_0119 = 3*IT_0118;
    const ccomplex_t IT_0120 = 0.166666666666667*IT_0119;
    const ccomplex_t IT_0121 = IT_0108*IT_0120;
    const ccomplex_t IT_0122 = IT_0100*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = 0.166666666666667*IT_0013 + -IT_0050 + -IT_0074
       + -IT_0099 + -IT_0123;
    const ccomplex_t IT_0125 = pow(m_sb_2, 2);
    const ccomplex_t IT_0126 = cpow((-2)*s_23 + IT_0006 + IT_0125 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0127 = V_tb*e_em*conj(V_Wp1)*U_st_01;
    const ccomplex_t IT_0128 = IT_0018*IT_0127;
    const ccomplex_t IT_0129 = m_t*V_tb*conj(V_u1)*e_em*IT_0020*U_st_11;
    const ccomplex_t IT_0130 = IT_0040*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0128 + (-0.5)*IT_0131);
    const ccomplex_t IT_0133 = -IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0018*IT_0020*IT_0039*U_sb_01;
    const ccomplex_t IT_0135 = 0.5*IT_0134;
    const ccomplex_t IT_0136 = IT_0133*IT_0135;
    const ccomplex_t IT_0137 = IT_0126*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = cpow((-2)*s_23 + IT_0006 + IT_0125 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0140 = V_tb*e_em*conj(V_Wp2)*U_st_01;
    const ccomplex_t IT_0141 = IT_0018*IT_0140;
    const ccomplex_t IT_0142 = m_t*V_tb*conj(V_u2)*e_em*IT_0020*U_st_11;
    const ccomplex_t IT_0143 = IT_0040*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0141 + (-0.5)*IT_0144);
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0018*IT_0020*IT_0039*U_sb_01;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = IT_0146*IT_0148;
    const ccomplex_t IT_0150 = IT_0139*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = -IT_0138 + -IT_0151;
    const ccomplex_t IT_0153 = (-0.5)*IT_0013;
    const ccomplex_t IT_0154 = IT_0152 + IT_0153;
    const ccomplex_t IT_0155 = conj(IT_0152) + conj(IT_0153);
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0158 = IT_0156*IT_0157;
    const ccomplex_t IT_0159 = IT_0008*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0162 = IT_0028*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0165 = IT_0018*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = m_t*conj(N_u1)*e_em*IT_0020*U_st_11;
    const ccomplex_t IT_0168 = IT_0040*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*(IT_0163 + 3*IT_0166 + 3
      *IT_0169);
    const ccomplex_t IT_0171 = 0.166666666666667*IT_0170;
    const ccomplex_t IT_0172 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0173 = IT_0028*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = m_b*N_d1*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0176 = IT_0019*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + 1.5*IT_0177);
    const ccomplex_t IT_0179 = (-0.333333333333333)*IT_0178;
    const ccomplex_t IT_0180 = IT_0171*IT_0179;
    const ccomplex_t IT_0181 = IT_0051*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0184 = IT_0028*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = m_b*N_d2*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0187 = IT_0019*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(IT_0185 + 1.5*IT_0188);
    const ccomplex_t IT_0190 = 0.666666666666667*IT_0189;
    const ccomplex_t IT_0191 = 1.5*IT_0190;
    const ccomplex_t IT_0192 = (-0.333333333333333)*IT_0191;
    const ccomplex_t IT_0193 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0194 = IT_0028*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0197 = IT_0018*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = m_t*conj(N_u2)*e_em*IT_0020*U_st_11;
    const ccomplex_t IT_0200 = IT_0040*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*(IT_0195 + 3*IT_0198 + 3
      *IT_0201);
    const ccomplex_t IT_0203 = 0.166666666666667*IT_0202;
    const ccomplex_t IT_0204 = IT_0192*IT_0203;
    const ccomplex_t IT_0205 = IT_0075*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0208 = IT_0028*IT_0207;
    const ccomplex_t IT_0209 = 1.4142135623731*IT_0208;
    const ccomplex_t IT_0210 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0211 = IT_0018*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = m_t*conj(N_u3)*e_em*IT_0020*U_st_11;
    const ccomplex_t IT_0214 = IT_0040*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*(IT_0209 + 3*IT_0212 + 3
      *IT_0215);
    const ccomplex_t IT_0217 = 0.166666666666667*IT_0216;
    const ccomplex_t IT_0218 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0219 = IT_0028*IT_0218;
    const ccomplex_t IT_0220 = 1.4142135623731*IT_0219;
    const ccomplex_t IT_0221 = m_b*N_d3*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0222 = IT_0019*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*(IT_0220 + 1.5*IT_0223);
    const ccomplex_t IT_0225 = 0.666666666666667*IT_0224;
    const ccomplex_t IT_0226 = 1.5*IT_0225;
    const ccomplex_t IT_0227 = 0.666666666666667*IT_0226;
    const ccomplex_t IT_0228 = 1.5*IT_0227;
    const ccomplex_t IT_0229 = 0.666666666666667*IT_0228;
    const ccomplex_t IT_0230 = 1.5*IT_0229;
    const ccomplex_t IT_0231 = (-0.333333333333333)*IT_0230;
    const ccomplex_t IT_0232 = IT_0217*IT_0231;
    const ccomplex_t IT_0233 = IT_0100*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0236 = IT_0028*IT_0235;
    const ccomplex_t IT_0237 = 1.4142135623731*IT_0236;
    const ccomplex_t IT_0238 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0239 = IT_0018*IT_0238;
    const ccomplex_t IT_0240 = 1.4142135623731*IT_0239;
    const ccomplex_t IT_0241 = m_t*conj(N_u4)*e_em*IT_0020*U_st_11;
    const ccomplex_t IT_0242 = IT_0040*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*(IT_0237 + 3*IT_0240 + 3
      *IT_0243);
    const ccomplex_t IT_0245 = 0.166666666666667*IT_0244;
    const ccomplex_t IT_0246 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0247 = IT_0028*IT_0246;
    const ccomplex_t IT_0248 = 1.4142135623731*IT_0247;
    const ccomplex_t IT_0249 = m_b*N_d4*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0250 = IT_0019*IT_0249;
    const ccomplex_t IT_0251 = 1.4142135623731*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*(IT_0248 + 1.5*IT_0251);
    const ccomplex_t IT_0253 = 0.666666666666667*IT_0252;
    const ccomplex_t IT_0254 = 1.5*IT_0253;
    const ccomplex_t IT_0255 = (-0.333333333333333)*IT_0254;
    const ccomplex_t IT_0256 = IT_0245*IT_0255;
    const ccomplex_t IT_0257 = IT_0014*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = 0.166666666666667*IT_0160 + -IT_0182 + -IT_0206
       + -IT_0234 + -IT_0258;
    const ccomplex_t IT_0260 = (-0.5)*IT_0160;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0016*IT_0018*IT_0020*U_st_01;
    const ccomplex_t IT_0262 = 0.5*IT_0261;
    const ccomplex_t IT_0263 = m_b*conj(U_d2)*V_tb*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0264 = IT_0019*IT_0263;
    const ccomplex_t IT_0265 = 1.4142135623731*IT_0264;
    const ccomplex_t IT_0266 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0267 = IT_0018*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*(IT_0265 + (-2)*IT_0267);
    const ccomplex_t IT_0269 = (-0.5)*IT_0268;
    const ccomplex_t IT_0270 = IT_0262*IT_0269;
    const ccomplex_t IT_0271 = IT_0139*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0016*IT_0018*IT_0020*U_st_01;
    const ccomplex_t IT_0274 = 0.5*IT_0273;
    const ccomplex_t IT_0275 = m_b*conj(U_d1)*V_tb*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0276 = IT_0019*IT_0275;
    const ccomplex_t IT_0277 = 1.4142135623731*IT_0276;
    const ccomplex_t IT_0278 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0279 = IT_0018*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*(IT_0277 + (-2)*IT_0279);
    const ccomplex_t IT_0281 = (-0.5)*IT_0280;
    const ccomplex_t IT_0282 = IT_0274*IT_0281;
    const ccomplex_t IT_0283 = IT_0126*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = IT_0272 + IT_0284;
    const ccomplex_t IT_0286 = IT_0260 + IT_0285;
    const ccomplex_t IT_0287 = conj(IT_0260) + conj(IT_0285);
    const ccomplex_t IT_0288 = (-18)*IT_0001;
    const ccomplex_t IT_0289 = 36*IT_0003;
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = m_b*m_t;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0293 = IT_0071*IT_0171;
    const ccomplex_t IT_0294 = IT_0051*IT_0292*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0296 = IT_0120*IT_0217;
    const ccomplex_t IT_0297 = IT_0100*IT_0295*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0299 = IT_0034*IT_0245;
    const ccomplex_t IT_0300 = IT_0014*IT_0298*IT_0299;
    const ccomplex_t IT_0301 = -IT_0294 + -IT_0297 + -IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0303 = IT_0148*IT_0262;
    const ccomplex_t IT_0304 = IT_0139*IT_0302*IT_0303;
    const ccomplex_t IT_0305 = m_t*IT_0160;
    const ccomplex_t IT_0306 = IT_0304 + (-0.5)*IT_0305;
    const ccomplex_t IT_0307 = 6*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0309 = IT_0096*IT_0192;
    const ccomplex_t IT_0310 = IT_0075*IT_0308*IT_0309;
    const ccomplex_t IT_0311 = m_t*IT_0182;
    const ccomplex_t IT_0312 = m_t*IT_0206;
    const ccomplex_t IT_0313 = m_t*IT_0234;
    const ccomplex_t IT_0314 = m_t*IT_0258;
    const ccomplex_t IT_0315 = 0.166666666666667*IT_0305 + IT_0310 + -IT_0311 
      + -IT_0312 + -IT_0313 + -IT_0314;
    const ccomplex_t IT_0316 = 18*IT_0315;
    const ccomplex_t IT_0317 = IT_0009*IT_0157;
    const ccomplex_t IT_0318 = IT_0008*IT_0317;
    const ccomplex_t IT_0319 = (0 + _Complex_I*1)*IT_0318;
    const ccomplex_t IT_0320 = m_sG*IT_0319;
    const ccomplex_t IT_0321 = IT_0047*IT_0255;
    const ccomplex_t IT_0322 = IT_0014*IT_0298*IT_0321;
    const ccomplex_t IT_0323 = IT_0059*IT_0179;
    const ccomplex_t IT_0324 = IT_0051*IT_0292*IT_0323;
    const ccomplex_t IT_0325 = IT_0108*IT_0231;
    const ccomplex_t IT_0326 = IT_0100*IT_0295*IT_0325;
    const ccomplex_t IT_0327 = 0.166666666666667*IT_0320 + -IT_0322 + -IT_0324
       + -IT_0326;
    const ccomplex_t IT_0328 = (-18)*IT_0327;
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0330 = IT_0135*IT_0274;
    const ccomplex_t IT_0331 = IT_0126*IT_0329*IT_0330;
    const ccomplex_t IT_0332 = m_b*IT_0138;
    const ccomplex_t IT_0333 = m_b*IT_0151;
    const ccomplex_t IT_0334 = (-0.5)*IT_0320 + -IT_0331 + -IT_0332 + -IT_0333;
    const ccomplex_t IT_0335 = (-6)*IT_0334;
    const ccomplex_t IT_0336 = IT_0307 + IT_0316 + IT_0328 + IT_0335;
    const ccomplex_t IT_0337 = 6*conj(IT_0306);
    const ccomplex_t IT_0338 = 18*conj(IT_0315);
    const ccomplex_t IT_0339 = (-18)*conj(IT_0327);
    const ccomplex_t IT_0340 = (-6)*conj(IT_0334);
    const ccomplex_t IT_0341 = IT_0337 + IT_0338 + IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = m_b*IT_0272;
    const ccomplex_t IT_0343 = IT_0146*IT_0269;
    const ccomplex_t IT_0344 = IT_0139*IT_0302*IT_0343;
    const ccomplex_t IT_0345 = m_b*IT_0284;
    const ccomplex_t IT_0346 = IT_0342 + IT_0344 + IT_0345;
    const ccomplex_t IT_0347 = 6*IT_0315;
    const ccomplex_t IT_0348 = 18*IT_0306;
    const ccomplex_t IT_0349 = (-18)*IT_0334;
    const ccomplex_t IT_0350 = (-6)*IT_0327;
    const ccomplex_t IT_0351 = IT_0347 + IT_0348 + IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = 6*conj(IT_0315);
    const ccomplex_t IT_0353 = 18*conj(IT_0306);
    const ccomplex_t IT_0354 = (-18)*conj(IT_0334);
    const ccomplex_t IT_0355 = (-6)*conj(IT_0327);
    const ccomplex_t IT_0356 = IT_0352 + IT_0353 + IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = m_t*IT_0013;
    const ccomplex_t IT_0358 = IT_0010*IT_0156;
    const ccomplex_t IT_0359 = IT_0008*IT_0358;
    const ccomplex_t IT_0360 = (0 + _Complex_I*1)*IT_0359;
    const ccomplex_t IT_0361 = m_sG*IT_0360;
    const ccomplex_t IT_0362 = IT_0133*IT_0281;
    const ccomplex_t IT_0363 = IT_0126*IT_0329*IT_0362;
    const ccomplex_t IT_0364 = (-0.5)*IT_0357 + 0.5*IT_0361 + -IT_0363;
    const ccomplex_t IT_0365 = 6*IT_0327;
    const ccomplex_t IT_0366 = 18*IT_0334;
    const ccomplex_t IT_0367 = (-18)*IT_0306;
    const ccomplex_t IT_0368 = (-6)*IT_0315;
    const ccomplex_t IT_0369 = IT_0365 + IT_0366 + IT_0367 + IT_0368;
    const ccomplex_t IT_0370 = 6*conj(IT_0327);
    const ccomplex_t IT_0371 = 18*conj(IT_0334);
    const ccomplex_t IT_0372 = (-18)*conj(IT_0306);
    const ccomplex_t IT_0373 = (-6)*conj(IT_0315);
    const ccomplex_t IT_0374 = IT_0370 + IT_0371 + IT_0372 + IT_0373;
    const ccomplex_t IT_0375 = IT_0088*IT_0203;
    const ccomplex_t IT_0376 = IT_0075*IT_0308*IT_0375;
    const ccomplex_t IT_0377 = m_t*IT_0074;
    const ccomplex_t IT_0378 = m_t*IT_0099;
    const ccomplex_t IT_0379 = m_t*IT_0123;
    const ccomplex_t IT_0380 = m_t*IT_0050;
    const ccomplex_t IT_0381 = 0.166666666666667*IT_0357 + (-0.166666666666667
      )*IT_0361 + IT_0376 + -IT_0377 + -IT_0378 + -IT_0379 + -IT_0380;
    const ccomplex_t IT_0382 = 6*IT_0334;
    const ccomplex_t IT_0383 = 18*IT_0327;
    const ccomplex_t IT_0384 = (-18)*IT_0315;
    const ccomplex_t IT_0385 = (-6)*IT_0306;
    const ccomplex_t IT_0386 = IT_0382 + IT_0383 + IT_0384 + IT_0385;
    const ccomplex_t IT_0387 = 6*conj(IT_0334);
    const ccomplex_t IT_0388 = 18*conj(IT_0327);
    const ccomplex_t IT_0389 = (-18)*conj(IT_0315);
    const ccomplex_t IT_0390 = (-6)*conj(IT_0306);
    const ccomplex_t IT_0391 = IT_0387 + IT_0388 + IT_0389 + IT_0390;
    const ccomplex_t IT_0392 = 6*s_34;
    const ccomplex_t IT_0393 = 18*s_34;
    const ccomplex_t IT_0394 = (-18)*s_34;
    const ccomplex_t IT_0395 = (-6)*s_34;
    const ccomplex_t IT_0396 = m_t*s_14;
    const ccomplex_t IT_0397 = 6*IT_0301;
    const ccomplex_t IT_0398 = 18*IT_0346;
    const ccomplex_t IT_0399 = (-18)*IT_0364;
    const ccomplex_t IT_0400 = (-6)*IT_0381;
    const ccomplex_t IT_0401 = 6*conj(IT_0301);
    const ccomplex_t IT_0402 = 18*conj(IT_0346);
    const ccomplex_t IT_0403 = (-18)*conj(IT_0364);
    const ccomplex_t IT_0404 = (-6)*conj(IT_0381);
    const ccomplex_t IT_0405 = IT_0287*IT_0369 + IT_0286*IT_0374 + conj
      (IT_0124)*(18*IT_0301 + 6*IT_0346 + (-6)*IT_0364 + (-18)*IT_0381) +
       IT_0124*(18*conj(IT_0301) + 6*conj(IT_0346) + (-6)*conj(IT_0364) + (-18)
      *conj(IT_0381)) + conj(IT_0259)*IT_0386 + IT_0259*IT_0391 + IT_0155*
      (IT_0397 + IT_0398 + IT_0399 + IT_0400) + IT_0154*(IT_0401 + IT_0402 +
       IT_0403 + IT_0404);
    const ccomplex_t IT_0406 = m_b*s_13;
    const ccomplex_t IT_0407 = 6*IT_0381;
    const ccomplex_t IT_0408 = 18*IT_0364;
    const ccomplex_t IT_0409 = (-18)*IT_0346;
    const ccomplex_t IT_0410 = (-6)*IT_0301;
    const ccomplex_t IT_0411 = 6*conj(IT_0381);
    const ccomplex_t IT_0412 = 18*conj(IT_0364);
    const ccomplex_t IT_0413 = (-18)*conj(IT_0346);
    const ccomplex_t IT_0414 = (-6)*conj(IT_0301);
    const ccomplex_t IT_0415 = conj(IT_0124)*IT_0336 + IT_0124*IT_0341 +
       IT_0155*IT_0351 + IT_0154*IT_0356 + conj(IT_0259)*((-18)*IT_0301 + (-6)
      *IT_0346 + 6*IT_0364 + 18*IT_0381) + IT_0259*((-18)*conj(IT_0301) + (-6)
      *conj(IT_0346) + 6*conj(IT_0364) + 18*conj(IT_0381)) + IT_0287*(IT_0407 +
       IT_0408 + IT_0409 + IT_0410) + IT_0286*(IT_0411 + IT_0412 + IT_0413 +
       IT_0414);
    const ccomplex_t IT_0416 = m_b*m_t*IT_0000;
    const ccomplex_t IT_0417 = (-18)*IT_0153;
    const ccomplex_t IT_0418 = (-6)*IT_0124;
    const ccomplex_t IT_0419 = (-18)*conj(IT_0153);
    const ccomplex_t IT_0420 = (-6)*conj(IT_0124);
    const ccomplex_t IT_0421 = ((-18)*conj(IT_0124) + (-6)*conj(IT_0153))
      *IT_0259 + ((-18)*IT_0124 + (-6)*IT_0153)*conj(IT_0259) + (-6)*conj
      (IT_0152)*(IT_0259 + 3*IT_0260 + 3*IT_0285) + (-6)*IT_0152*(conj(IT_0259) 
      + 3*conj(IT_0260) + 3*conj(IT_0285)) + IT_0287*(IT_0417 + IT_0418) +
       IT_0286*(IT_0419 + IT_0420);
    const ccomplex_t IT_0422 = IT_0005*(conj(IT_0124)*IT_0154 + IT_0124
      *IT_0155 + conj(IT_0259)*IT_0286 + IT_0259*IT_0287) + (IT_0124*conj
      (IT_0124) + IT_0154*IT_0155 + IT_0259*conj(IT_0259) + IT_0286*IT_0287)
      *IT_0290 + IT_0291*(conj(IT_0301)*IT_0336 + IT_0301*IT_0341 + conj(IT_0346
      )*IT_0351 + IT_0346*IT_0356 + conj(IT_0364)*IT_0369 + IT_0364*IT_0374 +
       conj(IT_0381)*IT_0386 + IT_0381*IT_0391) + (conj(IT_0306)*IT_0315 +
       IT_0306*conj(IT_0315) + conj(IT_0327)*IT_0334 + IT_0327*conj(IT_0334) +
       conj(IT_0301)*IT_0346 + IT_0301*conj(IT_0346) + conj(IT_0364)*IT_0381 +
       IT_0364*conj(IT_0381))*IT_0392 + (IT_0301*conj(IT_0301) + IT_0306*conj
      (IT_0306) + IT_0315*conj(IT_0315) + IT_0327*conj(IT_0327) + IT_0334*conj
      (IT_0334) + IT_0346*conj(IT_0346) + IT_0364*conj(IT_0364) + IT_0381*conj
      (IT_0381))*IT_0393 + (conj(IT_0315)*IT_0327 + IT_0315*conj(IT_0327) + conj
      (IT_0306)*IT_0334 + IT_0306*conj(IT_0334) + conj(IT_0346)*IT_0364 +
       IT_0346*conj(IT_0364) + conj(IT_0301)*IT_0381 + IT_0301*conj(IT_0381))
      *IT_0394 + (conj(IT_0306)*IT_0327 + IT_0306*conj(IT_0327) + conj(IT_0315)
      *IT_0334 + IT_0315*conj(IT_0334) + conj(IT_0301)*IT_0364 + IT_0301*conj
      (IT_0364) + conj(IT_0346)*IT_0381 + IT_0346*conj(IT_0381))*IT_0395 +
       IT_0396*IT_0405 + IT_0406*IT_0415 + IT_0416*IT_0421;
    return create_ccomplex_return(IT_0422);
}

