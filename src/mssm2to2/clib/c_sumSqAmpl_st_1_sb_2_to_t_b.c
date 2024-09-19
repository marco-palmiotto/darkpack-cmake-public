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
    const ccomplex_t IT_0068 = 3*IT_0067;
    const ccomplex_t IT_0069 = 0.166666666666667*IT_0068;
    const ccomplex_t IT_0070 = IT_0057*IT_0069;
    const ccomplex_t IT_0071 = IT_0049*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0074 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0075 = IT_0016*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = m_t*N_u3*e_em*IT_0025*U_st_00;
    const ccomplex_t IT_0078 = IT_0024*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0076 + (-0.75)*IT_0079);
    const ccomplex_t IT_0081 = 0.666666666666667*IT_0080;
    const ccomplex_t IT_0082 = m_b*conj(N_d3)*e_em*IT_0025*U_sb_11;
    const ccomplex_t IT_0083 = IT_0033*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0086 = IT_0016*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0089 = IT_0023*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0084 + 0.333333333333333
      *IT_0087 + -IT_0090);
    const ccomplex_t IT_0092 = 3*IT_0091;
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = IT_0081*IT_0093;
    const ccomplex_t IT_0095 = IT_0073*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0098 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0099 = IT_0016*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = m_t*N_u4*e_em*IT_0025*U_st_00;
    const ccomplex_t IT_0102 = IT_0024*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + (-0.75)*IT_0103);
    const ccomplex_t IT_0105 = 0.666666666666667*IT_0104;
    const ccomplex_t IT_0106 = m_b*conj(N_d4)*e_em*IT_0025*U_sb_11;
    const ccomplex_t IT_0107 = IT_0033*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0110 = IT_0023*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0113 = IT_0016*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0108 + -IT_0111 +
       0.333333333333333*IT_0114);
    const ccomplex_t IT_0116 = 3*IT_0115;
    const ccomplex_t IT_0117 = 0.166666666666667*IT_0116;
    const ccomplex_t IT_0118 = IT_0105*IT_0117;
    const ccomplex_t IT_0119 = IT_0097*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = 0.166666666666667*IT_0013 + -IT_0048 + -IT_0072
       + -IT_0096 + -IT_0120;
    const ccomplex_t IT_0122 = (-0.5)*IT_0013;
    const ccomplex_t IT_0123 = pow(m_sb_2, 2);
    const ccomplex_t IT_0124 = cpow((-2)*s_23 + IT_0006 + IT_0123 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0125 = V_tb*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0126 = IT_0023*IT_0125;
    const ccomplex_t IT_0127 = m_t*V_tb*conj(V_u1)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0128 = IT_0024*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + (-0.5)*IT_0129);
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0021*IT_0023*IT_0025*U_sb_01;
    const ccomplex_t IT_0133 = 0.5*IT_0132;
    const ccomplex_t IT_0134 = IT_0131*IT_0133;
    const ccomplex_t IT_0135 = IT_0124*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = cpow((-2)*s_23 + IT_0006 + IT_0123 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0138 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0139 = IT_0023*IT_0138;
    const ccomplex_t IT_0140 = m_t*V_tb*conj(V_u2)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0141 = IT_0024*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + (-0.5)*IT_0142);
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0021*IT_0023*IT_0025*U_sb_01;
    const ccomplex_t IT_0146 = 0.5*IT_0145;
    const ccomplex_t IT_0147 = IT_0144*IT_0146;
    const ccomplex_t IT_0148 = IT_0137*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = -IT_0136 + -IT_0149;
    const ccomplex_t IT_0151 = IT_0122 + IT_0150;
    const ccomplex_t IT_0152 = conj(IT_0122) + conj(IT_0150);
    const ccomplex_t IT_0153 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0155 = IT_0153*IT_0154;
    const ccomplex_t IT_0156 = IT_0008*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0159 = IT_0016*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0162 = IT_0023*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = m_t*conj(N_u1)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0165 = IT_0024*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0160 + 3*IT_0163 + 3
      *IT_0166);
    const ccomplex_t IT_0168 = 0.166666666666667*IT_0167;
    const ccomplex_t IT_0169 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0170 = IT_0016*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = m_b*N_d1*e_em*IT_0025*U_sb_01;
    const ccomplex_t IT_0173 = IT_0033*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0171 + 1.5*IT_0174);
    const ccomplex_t IT_0176 = (-0.333333333333333)*IT_0175;
    const ccomplex_t IT_0177 = IT_0168*IT_0176;
    const ccomplex_t IT_0178 = IT_0014*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0181 = IT_0016*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0184 = IT_0023*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = m_t*conj(N_u2)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0187 = IT_0024*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(IT_0182 + 3*IT_0185 + 3
      *IT_0188);
    const ccomplex_t IT_0190 = 0.166666666666667*IT_0189;
    const ccomplex_t IT_0191 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0192 = IT_0016*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = m_b*N_d2*e_em*IT_0025*U_sb_01;
    const ccomplex_t IT_0195 = IT_0033*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + 1.5*IT_0196);
    const ccomplex_t IT_0198 = (-0.333333333333333)*IT_0197;
    const ccomplex_t IT_0199 = IT_0190*IT_0198;
    const ccomplex_t IT_0200 = IT_0049*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0203 = IT_0016*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0206 = IT_0023*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = m_t*conj(N_u3)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0209 = IT_0024*IT_0208;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*(IT_0204 + 3*IT_0207 + 3
      *IT_0210);
    const ccomplex_t IT_0212 = 0.166666666666667*IT_0211;
    const ccomplex_t IT_0213 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0214 = IT_0016*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = m_b*N_d3*e_em*IT_0025*U_sb_01;
    const ccomplex_t IT_0217 = IT_0033*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*(IT_0215 + 1.5*IT_0218);
    const ccomplex_t IT_0220 = (-0.333333333333333)*IT_0219;
    const ccomplex_t IT_0221 = IT_0212*IT_0220;
    const ccomplex_t IT_0222 = IT_0073*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0225 = IT_0016*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0228 = IT_0023*IT_0227;
    const ccomplex_t IT_0229 = 1.4142135623731*IT_0228;
    const ccomplex_t IT_0230 = m_t*conj(N_u4)*e_em*IT_0025*U_st_10;
    const ccomplex_t IT_0231 = IT_0024*IT_0230;
    const ccomplex_t IT_0232 = 1.4142135623731*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*(IT_0226 + 3*IT_0229 + 3
      *IT_0232);
    const ccomplex_t IT_0234 = 0.166666666666667*IT_0233;
    const ccomplex_t IT_0235 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0236 = IT_0016*IT_0235;
    const ccomplex_t IT_0237 = 1.4142135623731*IT_0236;
    const ccomplex_t IT_0238 = m_b*N_d4*e_em*IT_0025*U_sb_01;
    const ccomplex_t IT_0239 = IT_0033*IT_0238;
    const ccomplex_t IT_0240 = 1.4142135623731*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*(IT_0237 + 1.5*IT_0240);
    const ccomplex_t IT_0242 = (-0.333333333333333)*IT_0241;
    const ccomplex_t IT_0243 = IT_0234*IT_0242;
    const ccomplex_t IT_0244 = IT_0097*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*IT_0244;
    const ccomplex_t IT_0246 = 0.166666666666667*IT_0157 + -IT_0179 + -IT_0201
       + -IT_0223 + -IT_0245;
    const ccomplex_t IT_0247 = (-0.5)*IT_0157;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0023*IT_0025*IT_0032*U_st_00;
    const ccomplex_t IT_0249 = 0.5*IT_0248;
    const ccomplex_t IT_0250 = m_b*conj(U_d1)*V_tb*e_em*IT_0025*U_sb_11;
    const ccomplex_t IT_0251 = IT_0033*IT_0250;
    const ccomplex_t IT_0252 = 1.4142135623731*IT_0251;
    const ccomplex_t IT_0253 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0254 = IT_0023*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*(IT_0252 + (-2)*IT_0254);
    const ccomplex_t IT_0256 = (-0.5)*IT_0255;
    const ccomplex_t IT_0257 = IT_0249*IT_0256;
    const ccomplex_t IT_0258 = IT_0124*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0023*IT_0025*IT_0032*U_st_00;
    const ccomplex_t IT_0261 = 0.5*IT_0260;
    const ccomplex_t IT_0262 = m_b*conj(U_d2)*V_tb*e_em*IT_0025*U_sb_11;
    const ccomplex_t IT_0263 = IT_0033*IT_0262;
    const ccomplex_t IT_0264 = 1.4142135623731*IT_0263;
    const ccomplex_t IT_0265 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0266 = IT_0023*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*(IT_0264 + (-2)*IT_0266);
    const ccomplex_t IT_0268 = (-0.5)*IT_0267;
    const ccomplex_t IT_0269 = IT_0261*IT_0268;
    const ccomplex_t IT_0270 = IT_0137*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = IT_0259 + IT_0271;
    const ccomplex_t IT_0273 = IT_0247 + IT_0272;
    const ccomplex_t IT_0274 = conj(IT_0247) + conj(IT_0272);
    const ccomplex_t IT_0275 = (-18)*IT_0001;
    const ccomplex_t IT_0276 = 36*IT_0003;
    const ccomplex_t IT_0277 = IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = m_b*m_t;
    const ccomplex_t IT_0279 = IT_0010*IT_0153;
    const ccomplex_t IT_0280 = IT_0008*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = m_sG*IT_0281;
    const ccomplex_t IT_0283 = (-0.5)*IT_0282;
    const ccomplex_t IT_0284 = m_b*IT_0272;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0286 = IT_0131*IT_0256;
    const ccomplex_t IT_0287 = IT_0124*IT_0285*IT_0286;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0289 = IT_0030*IT_0176;
    const ccomplex_t IT_0290 = IT_0014*IT_0288*IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0292 = IT_0105*IT_0242;
    const ccomplex_t IT_0293 = IT_0097*IT_0291*IT_0292;
    const ccomplex_t IT_0294 = IT_0009*IT_0154;
    const ccomplex_t IT_0295 = IT_0008*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = m_sG*IT_0296;
    const ccomplex_t IT_0298 = m_t*IT_0157;
    const ccomplex_t IT_0299 = m_t*IT_0179;
    const ccomplex_t IT_0300 = m_t*IT_0201;
    const ccomplex_t IT_0301 = m_t*IT_0223;
    const ccomplex_t IT_0302 = m_t*IT_0245;
    const ccomplex_t IT_0303 = IT_0290 + IT_0293 + (-0.166666666666667)
      *IT_0297 + 0.166666666666667*IT_0298 + -IT_0299 + -IT_0300 + -IT_0301 + 
      -IT_0302;
    const ccomplex_t IT_0304 = 6*IT_0303;
    const ccomplex_t IT_0305 = IT_0133*IT_0249;
    const ccomplex_t IT_0306 = IT_0124*IT_0285*IT_0305;
    const ccomplex_t IT_0307 = (-0.5)*IT_0298 + IT_0306;
    const ccomplex_t IT_0308 = 18*IT_0307;
    const ccomplex_t IT_0309 = 0.5*IT_0297;
    const ccomplex_t IT_0310 = 18*IT_0309;
    const ccomplex_t IT_0311 = m_b*IT_0149;
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = (-18)*IT_0312;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0315 = IT_0146*IT_0261;
    const ccomplex_t IT_0316 = IT_0137*IT_0314*IT_0315;
    const ccomplex_t IT_0317 = m_b*IT_0136;
    const ccomplex_t IT_0318 = -IT_0316 + -IT_0317;
    const ccomplex_t IT_0319 = (-18)*IT_0318;
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0321 = IT_0081*IT_0220;
    const ccomplex_t IT_0322 = IT_0073*IT_0320*IT_0321;
    const ccomplex_t IT_0323 = -IT_0322;
    const ccomplex_t IT_0324 = (-6)*IT_0323;
    const ccomplex_t IT_0325 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0326 = IT_0057*IT_0198;
    const ccomplex_t IT_0327 = IT_0049*IT_0325*IT_0326;
    const ccomplex_t IT_0328 = -IT_0327;
    const ccomplex_t IT_0329 = (-6)*IT_0328;
    const ccomplex_t IT_0330 = IT_0304 + IT_0308 + IT_0310 + IT_0313 + IT_0319
       + IT_0324 + IT_0329;
    const ccomplex_t IT_0331 = 6*conj(IT_0303);
    const ccomplex_t IT_0332 = 18*conj(IT_0307);
    const ccomplex_t IT_0333 = 18*conj(IT_0309);
    const ccomplex_t IT_0334 = (-18)*conj(IT_0312);
    const ccomplex_t IT_0335 = (-18)*conj(IT_0318);
    const ccomplex_t IT_0336 = (-6)*conj(IT_0323);
    const ccomplex_t IT_0337 = (-6)*conj(IT_0328);
    const ccomplex_t IT_0338 = IT_0331 + IT_0332 + IT_0333 + IT_0334 + IT_0335
       + IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = IT_0045*IT_0168;
    const ccomplex_t IT_0340 = IT_0014*IT_0288*IT_0339;
    const ccomplex_t IT_0341 = IT_0117*IT_0234;
    const ccomplex_t IT_0342 = IT_0097*IT_0291*IT_0341;
    const ccomplex_t IT_0343 = m_t*IT_0013;
    const ccomplex_t IT_0344 = m_t*IT_0048;
    const ccomplex_t IT_0345 = m_t*IT_0072;
    const ccomplex_t IT_0346 = m_t*IT_0096;
    const ccomplex_t IT_0347 = m_t*IT_0120;
    const ccomplex_t IT_0348 = IT_0340 + IT_0342 + 0.166666666666667*IT_0343 +
       -IT_0344 + -IT_0345 + -IT_0346 + -IT_0347;
    const ccomplex_t IT_0349 = 6*IT_0312;
    const ccomplex_t IT_0350 = 6*IT_0318;
    const ccomplex_t IT_0351 = 18*IT_0323;
    const ccomplex_t IT_0352 = 18*IT_0328;
    const ccomplex_t IT_0353 = (-18)*IT_0303;
    const ccomplex_t IT_0354 = (-6)*IT_0307;
    const ccomplex_t IT_0355 = (-6)*IT_0309;
    const ccomplex_t IT_0356 = IT_0349 + IT_0350 + IT_0351 + IT_0352 + IT_0353
       + IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = 6*conj(IT_0312);
    const ccomplex_t IT_0358 = 6*conj(IT_0318);
    const ccomplex_t IT_0359 = 18*conj(IT_0323);
    const ccomplex_t IT_0360 = 18*conj(IT_0328);
    const ccomplex_t IT_0361 = (-18)*conj(IT_0303);
    const ccomplex_t IT_0362 = (-6)*conj(IT_0307);
    const ccomplex_t IT_0363 = (-6)*conj(IT_0309);
    const ccomplex_t IT_0364 = IT_0357 + IT_0358 + IT_0359 + IT_0360 + IT_0361
       + IT_0362 + IT_0363;
    const ccomplex_t IT_0365 = IT_0144*IT_0268;
    const ccomplex_t IT_0366 = IT_0137*IT_0314*IT_0365;
    const ccomplex_t IT_0367 = (-0.5)*IT_0343 + -IT_0366;
    const ccomplex_t IT_0368 = 6*IT_0323;
    const ccomplex_t IT_0369 = 6*IT_0328;
    const ccomplex_t IT_0370 = 18*IT_0312;
    const ccomplex_t IT_0371 = 18*IT_0318;
    const ccomplex_t IT_0372 = (-18)*IT_0307;
    const ccomplex_t IT_0373 = (-18)*IT_0309;
    const ccomplex_t IT_0374 = (-6)*IT_0303;
    const ccomplex_t IT_0375 = IT_0368 + IT_0369 + IT_0370 + IT_0371 + IT_0372
       + IT_0373 + IT_0374;
    const ccomplex_t IT_0376 = 6*conj(IT_0323);
    const ccomplex_t IT_0377 = 6*conj(IT_0328);
    const ccomplex_t IT_0378 = 18*conj(IT_0312);
    const ccomplex_t IT_0379 = 18*conj(IT_0318);
    const ccomplex_t IT_0380 = (-18)*conj(IT_0307);
    const ccomplex_t IT_0381 = (-18)*conj(IT_0309);
    const ccomplex_t IT_0382 = (-6)*conj(IT_0303);
    const ccomplex_t IT_0383 = IT_0376 + IT_0377 + IT_0378 + IT_0379 + IT_0380
       + IT_0381 + IT_0382;
    const ccomplex_t IT_0384 = IT_0093*IT_0212;
    const ccomplex_t IT_0385 = IT_0073*IT_0320*IT_0384;
    const ccomplex_t IT_0386 = -IT_0385;
    const ccomplex_t IT_0387 = 0.166666666666667*IT_0282;
    const ccomplex_t IT_0388 = IT_0069*IT_0190;
    const ccomplex_t IT_0389 = IT_0049*IT_0325*IT_0388;
    const ccomplex_t IT_0390 = -IT_0389;
    const ccomplex_t IT_0391 = 6*IT_0307;
    const ccomplex_t IT_0392 = 6*IT_0309;
    const ccomplex_t IT_0393 = 18*IT_0303;
    const ccomplex_t IT_0394 = (-18)*IT_0323;
    const ccomplex_t IT_0395 = (-18)*IT_0328;
    const ccomplex_t IT_0396 = (-6)*IT_0312;
    const ccomplex_t IT_0397 = (-6)*IT_0318;
    const ccomplex_t IT_0398 = IT_0391 + IT_0392 + IT_0393 + IT_0394 + IT_0395
       + IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = 6*conj(IT_0307);
    const ccomplex_t IT_0400 = 6*conj(IT_0309);
    const ccomplex_t IT_0401 = 18*conj(IT_0303);
    const ccomplex_t IT_0402 = (-18)*conj(IT_0323);
    const ccomplex_t IT_0403 = (-18)*conj(IT_0328);
    const ccomplex_t IT_0404 = (-6)*conj(IT_0312);
    const ccomplex_t IT_0405 = (-6)*conj(IT_0318);
    const ccomplex_t IT_0406 = IT_0399 + IT_0400 + IT_0401 + IT_0402 + IT_0403
       + IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = (conj(IT_0283) + conj(IT_0284) + conj(IT_0287))
      *IT_0330 + (IT_0283 + IT_0284 + IT_0287)*IT_0338 + conj(IT_0348)*IT_0356 +
       IT_0348*IT_0364 + conj(IT_0367)*IT_0375 + IT_0367*IT_0383 + (conj(IT_0386
      ) + conj(IT_0387) + conj(IT_0390))*IT_0398 + (IT_0386 + IT_0387 + IT_0390)
      *IT_0406;
    const ccomplex_t IT_0408 = conj(IT_0312) + conj(IT_0318);
    const ccomplex_t IT_0409 = IT_0323 + IT_0328;
    const ccomplex_t IT_0410 = conj(IT_0323) + conj(IT_0328);
    const ccomplex_t IT_0411 = conj(IT_0307) + conj(IT_0309);
    const ccomplex_t IT_0412 = 6*s_34;
    const ccomplex_t IT_0413 = IT_0312 + IT_0318;
    const ccomplex_t IT_0414 = conj(IT_0387) + conj(IT_0390);
    const ccomplex_t IT_0415 = conj(IT_0284) + conj(IT_0287);
    const ccomplex_t IT_0416 = 18*s_34;
    const ccomplex_t IT_0417 = IT_0312*IT_0412;
    const ccomplex_t IT_0418 = IT_0323*IT_0416;
    const ccomplex_t IT_0419 = IT_0307 + IT_0309;
    const ccomplex_t IT_0420 = (-18)*s_34;
    const ccomplex_t IT_0421 = conj(IT_0303)*IT_0420;
    const ccomplex_t IT_0422 = IT_0386*IT_0416;
    const ccomplex_t IT_0423 = IT_0348*IT_0420;
    const ccomplex_t IT_0424 = conj(IT_0348)*IT_0420;
    const ccomplex_t IT_0425 = (-6)*s_34;
    const ccomplex_t IT_0426 = conj(IT_0303)*IT_0425;
    const ccomplex_t IT_0427 = IT_0307*IT_0416;
    const ccomplex_t IT_0428 = IT_0312*IT_0420;
    const ccomplex_t IT_0429 = IT_0323*IT_0425;
    const ccomplex_t IT_0430 = conj(IT_0303)*IT_0412;
    const ccomplex_t IT_0431 = conj(IT_0312)*IT_0420;
    const ccomplex_t IT_0432 = conj(IT_0323)*IT_0425;
    const ccomplex_t IT_0433 = IT_0386*IT_0412;
    const ccomplex_t IT_0434 = IT_0387*IT_0412;
    const ccomplex_t IT_0435 = IT_0433 + IT_0434;
    const ccomplex_t IT_0436 = IT_0283*IT_0416;
    const ccomplex_t IT_0437 = IT_0367*IT_0420;
    const ccomplex_t IT_0438 = IT_0348*IT_0425;
    const ccomplex_t IT_0439 = IT_0390*IT_0412;
    const ccomplex_t IT_0440 = IT_0284*IT_0416;
    const ccomplex_t IT_0441 = conj(IT_0386)*IT_0412;
    const ccomplex_t IT_0442 = conj(IT_0387)*IT_0412;
    const ccomplex_t IT_0443 = IT_0441 + IT_0442;
    const ccomplex_t IT_0444 = conj(IT_0367)*IT_0420;
    const ccomplex_t IT_0445 = conj(IT_0348)*IT_0425;
    const ccomplex_t IT_0446 = IT_0284 + IT_0287;
    const ccomplex_t IT_0447 = conj(IT_0390)*IT_0412;
    const ccomplex_t IT_0448 = conj(IT_0283)*IT_0416;
    const ccomplex_t IT_0449 = IT_0387*IT_0416;
    const ccomplex_t IT_0450 = IT_0367*IT_0425;
    const ccomplex_t IT_0451 = conj(IT_0386)*IT_0416;
    const ccomplex_t IT_0452 = conj(IT_0367)*IT_0425;
    const ccomplex_t IT_0453 = m_t*s_14;
    const ccomplex_t IT_0454 = 6*IT_0121;
    const ccomplex_t IT_0455 = 18*IT_0122;
    const ccomplex_t IT_0456 = 18*IT_0150;
    const ccomplex_t IT_0457 = 6*conj(IT_0121);
    const ccomplex_t IT_0458 = 18*conj(IT_0122);
    const ccomplex_t IT_0459 = 18*conj(IT_0150);
    const ccomplex_t IT_0460 = 6*IT_0386;
    const ccomplex_t IT_0461 = 6*IT_0387;
    const ccomplex_t IT_0462 = 6*IT_0390;
    const ccomplex_t IT_0463 = 18*IT_0283;
    const ccomplex_t IT_0464 = (-18)*IT_0367;
    const ccomplex_t IT_0465 = (-6)*IT_0348;
    const ccomplex_t IT_0466 = 6*conj(IT_0386);
    const ccomplex_t IT_0467 = 6*conj(IT_0387);
    const ccomplex_t IT_0468 = 6*conj(IT_0390);
    const ccomplex_t IT_0469 = 18*conj(IT_0283);
    const ccomplex_t IT_0470 = (-18)*conj(IT_0367);
    const ccomplex_t IT_0471 = (-6)*conj(IT_0348);
    const ccomplex_t IT_0472 = conj(IT_0246)*IT_0356 + IT_0246*IT_0364 +
       IT_0274*IT_0375 + IT_0273*IT_0383 + 6*conj(IT_0121)*(IT_0283 + (-3)
      *IT_0348 + -IT_0367 + 3*IT_0386 + 3*IT_0387 + 3*IT_0390) + 6*IT_0121*(conj
      (IT_0283) + (-3)*conj(IT_0348) + -conj(IT_0367) + 3*conj(IT_0386) + 3*conj
      (IT_0387) + 3*conj(IT_0390)) + IT_0415*(IT_0454 + IT_0455 + IT_0456) +
       IT_0446*(IT_0457 + IT_0458 + IT_0459) + IT_0152*(IT_0460 + IT_0461 +
       IT_0462 + IT_0463 + IT_0464 + IT_0465) + IT_0151*(IT_0466 + IT_0467 +
       IT_0468 + IT_0469 + IT_0470 + IT_0471);
    const ccomplex_t IT_0473 = m_b*s_13;
    const ccomplex_t IT_0474 = (-18)*IT_0247;
    const ccomplex_t IT_0475 = (-18)*IT_0272;
    const ccomplex_t IT_0476 = (-6)*IT_0246;
    const ccomplex_t IT_0477 = IT_0474 + IT_0475 + IT_0476;
    const ccomplex_t IT_0478 = (-18)*conj(IT_0247);
    const ccomplex_t IT_0479 = (-18)*conj(IT_0272);
    const ccomplex_t IT_0480 = (-6)*conj(IT_0246);
    const ccomplex_t IT_0481 = IT_0478 + IT_0479 + IT_0480;
    const ccomplex_t IT_0482 = 6*IT_0348;
    const ccomplex_t IT_0483 = 18*IT_0367;
    const ccomplex_t IT_0484 = (-18)*IT_0283;
    const ccomplex_t IT_0485 = (-6)*IT_0386;
    const ccomplex_t IT_0486 = (-6)*IT_0387;
    const ccomplex_t IT_0487 = (-6)*IT_0390;
    const ccomplex_t IT_0488 = 6*conj(IT_0348);
    const ccomplex_t IT_0489 = 18*conj(IT_0367);
    const ccomplex_t IT_0490 = (-18)*conj(IT_0283);
    const ccomplex_t IT_0491 = (-6)*conj(IT_0386);
    const ccomplex_t IT_0492 = (-6)*conj(IT_0387);
    const ccomplex_t IT_0493 = (-6)*conj(IT_0390);
    const ccomplex_t IT_0494 = IT_0152*IT_0330 + IT_0151*IT_0338 + (-6)*conj
      (IT_0246)*(IT_0283 + (-3)*IT_0348 + -IT_0367 + 3*IT_0386 + 3*IT_0387 + 3
      *IT_0390) + (-6)*IT_0246*(conj(IT_0283) + (-3)*conj(IT_0348) + -conj
      (IT_0367) + 3*conj(IT_0386) + 3*conj(IT_0387) + 3*conj(IT_0390)) + conj
      (IT_0121)*IT_0398 + IT_0121*IT_0406 + IT_0415*IT_0477 + IT_0446*IT_0481 +
       IT_0274*(IT_0482 + IT_0483 + IT_0484 + IT_0485 + IT_0486 + IT_0487) +
       IT_0273*(IT_0488 + IT_0489 + IT_0490 + IT_0491 + IT_0492 + IT_0493);
    const ccomplex_t IT_0495 = m_b*m_t*IT_0000;
    const ccomplex_t IT_0496 = (-18)*IT_0122;
    const ccomplex_t IT_0497 = (-6)*IT_0121;
    const ccomplex_t IT_0498 = (-18)*conj(IT_0122);
    const ccomplex_t IT_0499 = (-6)*conj(IT_0121);
    const ccomplex_t IT_0500 = IT_0495*(((-18)*conj(IT_0121) + (-6)*conj
      (IT_0122))*IT_0246 + ((-18)*IT_0121 + (-6)*IT_0122)*conj(IT_0246) + conj
      (IT_0150)*IT_0477 + IT_0150*IT_0481 + IT_0274*(IT_0496 + IT_0497) +
       IT_0273*(IT_0498 + IT_0499));
    const ccomplex_t IT_0501 = IT_0005*(conj(IT_0121)*IT_0151 + IT_0121
      *IT_0152 + conj(IT_0246)*IT_0273 + IT_0246*IT_0274) + (IT_0121*conj
      (IT_0121) + IT_0151*IT_0152 + IT_0246*conj(IT_0246) + IT_0273*IT_0274)
      *IT_0277 + IT_0278*IT_0407 + (conj(IT_0348)*IT_0367 + IT_0348*conj(IT_0367
      ) + conj(IT_0283)*IT_0390 + IT_0283*conj(IT_0390) + IT_0408*IT_0409 +
       IT_0318*IT_0410 + IT_0303*IT_0411)*IT_0412 + (IT_0303*conj(IT_0303) +
       IT_0348*conj(IT_0348) + IT_0367*conj(IT_0367) + IT_0328*IT_0410 + IT_0309
      *IT_0411 + IT_0408*IT_0413 + IT_0390*IT_0414 + IT_0287*IT_0415)*IT_0416 +
       IT_0410*(IT_0417 + IT_0418) + (IT_0303*IT_0410 + IT_0318*IT_0411 + conj
      (IT_0318)*IT_0419)*IT_0420 + IT_0409*IT_0421 + conj(IT_0386)*(IT_0422 +
       IT_0423) + IT_0386*IT_0424 + (conj(IT_0367)*IT_0386 + IT_0367*conj
      (IT_0386) + IT_0303*IT_0408 + IT_0328*IT_0411 + conj(IT_0328)*IT_0419)
      *IT_0425 + IT_0413*IT_0426 + IT_0411*(IT_0427 + IT_0428 + IT_0429) +
       IT_0419*(IT_0430 + IT_0431 + IT_0432) + conj(IT_0283)*(IT_0435 + IT_0436 
      + IT_0437 + IT_0438) + IT_0415*(IT_0435 + IT_0436 + IT_0437 + IT_0438 +
       IT_0439 + IT_0440) + IT_0283*(IT_0443 + IT_0444 + IT_0445) + IT_0446*
      (IT_0443 + IT_0444 + IT_0445 + IT_0447 + IT_0448) + IT_0414*(IT_0422 +
       IT_0423 + IT_0449 + IT_0450) + (IT_0387 + IT_0390)*(IT_0424 + IT_0451 +
       IT_0452) + IT_0453*IT_0472 + IT_0473*IT_0494 + IT_0500;
    return create_ccomplex_return(IT_0501);
}

