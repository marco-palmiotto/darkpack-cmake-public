#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_sb_2_to_anti_b_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_sb_2_to_anti_b_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
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
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 2*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = IT_0013 + IT_0016;
    const ccomplex_t IT_0018 = IT_0007 + IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = N_B1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0022 = IT_0009*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = N_W1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0025 = IT_0003*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = cpow(IT_0001, -1);
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = m_b*N_d1*e_em*IT_0029*conj(U_sb_11);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0023 + (-3)*IT_0026 + 3
      *IT_0032);
    const ccomplex_t IT_0034 = 0.166666666666667*IT_0033;
    const ccomplex_t IT_0035 = IT_0020*IT_0034;
    const ccomplex_t IT_0036 = pow(m_b, 2);
    const ccomplex_t IT_0037 = pow(m_sb_2, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0036 + IT_0037 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0000*IT_0035*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0041 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0042 = IT_0004*IT_0041;
    const ccomplex_t IT_0043 = 2*IT_0042;
    const ccomplex_t IT_0044 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0045 = IT_0010*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0048 = IT_0004*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = IT_0046 + IT_0049;
    const ccomplex_t IT_0051 = IT_0043 + IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0055 = IT_0009*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0058 = IT_0003*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = m_b*N_d3*e_em*IT_0029*conj(U_sb_11);
    const ccomplex_t IT_0061 = IT_0028*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0056 + (-3)*IT_0059 + 3
      *IT_0062);
    const ccomplex_t IT_0064 = 0.166666666666667*IT_0063;
    const ccomplex_t IT_0065 = IT_0053*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_23 + IT_0036 + IT_0037 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0040*IT_0065*IT_0066;
    const ccomplex_t IT_0068 = pow(m_C_1, 2);
    const ccomplex_t IT_0069 = cpow((-2)*s_13 + IT_0036 + IT_0068 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0070 = sin(beta);
    const ccomplex_t IT_0071 = cpow(IT_0070, -1);
    const ccomplex_t IT_0072 = m_t*conj(U_sb_01)*(m_t*IT_0001*IT_0071*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0027*IT_0070*conj(U_sb_01)*U_st_00 + conj
      (U_sb_11)*(mu_h*U_st_00 + m_t*IT_0027*IT_0071*U_st_10));
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0070*conj
      (U_sb_11)*U_st_00 + (-1.4142135623731)*e_em*IT_0003*((-0.5)*IT_0029
      *IT_0072 + m_W*IT_0001*IT_0070*conj(U_sb_01)*U_st_00) + A_t*IT_0001*conj
      (U_sb_01)*U_st_10);
    const ccomplex_t IT_0074 = IT_0003*IT_0071;
    const ccomplex_t IT_0075 = m_t*V_tb*V_u1*e_em*IT_0029*conj(U_st_10);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0081 = IT_0003*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = IT_0079 + IT_0084;
    const ccomplex_t IT_0086 = IT_0073*IT_0085;
    const ccomplex_t IT_0087 = IT_0069*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_13 + IT_0036 + IT_0068 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0090 = m_t*conj(U_sb_01)*(m_t*IT_0001*IT_0071*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0027*IT_0070*conj(U_sb_01)*U_st_01 + conj
      (U_sb_11)*(mu_h*U_st_01 + m_t*IT_0027*IT_0071*U_st_11));
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0070*conj
      (U_sb_11)*U_st_01 + (-1.4142135623731)*e_em*IT_0003*((-0.5)*IT_0029
      *IT_0090 + m_W*IT_0001*IT_0070*conj(U_sb_01)*U_st_01) + A_t*IT_0001*conj
      (U_sb_01)*U_st_11);
    const ccomplex_t IT_0092 = m_t*V_tb*V_u1*e_em*IT_0029*conj(U_st_11);
    const ccomplex_t IT_0093 = IT_0074*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0098 = IT_0003*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = IT_0096 + IT_0101;
    const ccomplex_t IT_0103 = IT_0091*IT_0102;
    const ccomplex_t IT_0104 = IT_0089*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_23 + IT_0036 + IT_0037 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0107 = conj(N_B2)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0108 = IT_0009*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = m_b*conj(N_d2)*e_em*IT_0029*conj(U_sb_01);
    const ccomplex_t IT_0111 = IT_0028*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + 1.5*IT_0112);
    const ccomplex_t IT_0114 = (-0.333333333333333)*IT_0113;
    const ccomplex_t IT_0115 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0116 = IT_0010*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0119 = IT_0004*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = IT_0117 + IT_0120;
    const ccomplex_t IT_0122 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0123 = IT_0004*IT_0122;
    const ccomplex_t IT_0124 = 2*IT_0123;
    const ccomplex_t IT_0125 = IT_0121 + IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = IT_0114*IT_0127;
    const ccomplex_t IT_0129 = IT_0106*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = m_b*IT_0130;
    const ccomplex_t IT_0132 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0133 = IT_0004*IT_0132;
    const ccomplex_t IT_0134 = 2*IT_0133;
    const ccomplex_t IT_0135 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0136 = IT_0010*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0139 = IT_0004*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = IT_0137 + IT_0140;
    const ccomplex_t IT_0142 = IT_0134 + IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = 0.5*IT_0143;
    const ccomplex_t IT_0145 = conj(N_B4)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0146 = IT_0009*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = m_b*conj(N_d4)*e_em*IT_0029*conj(U_sb_01);
    const ccomplex_t IT_0149 = IT_0028*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + 1.5*IT_0150);
    const ccomplex_t IT_0152 = (-0.333333333333333)*IT_0151;
    const ccomplex_t IT_0153 = IT_0144*IT_0152;
    const ccomplex_t IT_0154 = cpow((-2)*s_23 + IT_0036 + IT_0037 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0155 = IT_0153*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = m_b*IT_0156;
    const ccomplex_t IT_0158 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0159 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0003*IT_0029*IT_0071*IT_0158*V_ub_mod;
    const ccomplex_t IT_0160 = 0.5*IT_0159;
    const ccomplex_t IT_0161 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0162 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0003*IT_0029*IT_0161*(mu_h*conj(U_sb_01) + m_b*IT_0027*IT_0071*conj
      (U_sb_11))*V_ub_mod;
    const ccomplex_t IT_0163 = 0.5*IT_0162;
    const ccomplex_t IT_0164 = IT_0160*IT_0163;
    const ccomplex_t IT_0165 = cpow((-2)*s_13 + IT_0036 + IT_0068 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0166 = IT_0164*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0169 = N_B2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0170 = IT_0009*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = N_W2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0173 = IT_0003*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = m_b*N_d2*e_em*IT_0029*conj(U_sb_11);
    const ccomplex_t IT_0176 = IT_0028*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0171 + (-3)*IT_0174 + 3
      *IT_0177);
    const ccomplex_t IT_0179 = 0.166666666666667*IT_0178;
    const ccomplex_t IT_0180 = IT_0127*IT_0179;
    const ccomplex_t IT_0181 = IT_0106*IT_0168*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0183 = N_B4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0184 = IT_0009*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = N_W4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0187 = IT_0003*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = m_b*N_d4*e_em*IT_0029*conj(U_sb_11);
    const ccomplex_t IT_0190 = IT_0028*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*(IT_0185 + (-3)*IT_0188 + 3
      *IT_0191);
    const ccomplex_t IT_0193 = 0.166666666666667*IT_0192;
    const ccomplex_t IT_0194 = IT_0144*IT_0193;
    const ccomplex_t IT_0195 = IT_0154*IT_0182*IT_0194;
    const ccomplex_t IT_0196 = cpow((-2)*s_13 + IT_0036 + IT_0068 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0003*IT_0158
      *V_ub_mod;
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = pow(m_u, 2);
    const ccomplex_t IT_0200 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0003
      *IT_0161*(m_W*IT_0001*IT_0070*conj(U_sb_01) + (-0.5)*IT_0029*(IT_0001
      *IT_0071*IT_0199*conj(U_sb_01) + m_b*(m_b*IT_0027*IT_0070*conj(U_sb_01) +
       mu_h*conj(U_sb_11))))*V_ub_mod;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = IT_0198*IT_0201;
    const ccomplex_t IT_0203 = IT_0196*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = cpow((-2)*s_13 + IT_0036 + IT_0068 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*V_cb*e_em*V_Wp1*IT_0003;
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = pow(m_c, 2);
    const ccomplex_t IT_0209 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0003*(m_W*IT_0001*IT_0070*conj(U_sb_01) + (-0.5)*IT_0029*(IT_0001
      *IT_0071*IT_0208*conj(U_sb_01) + m_b*(m_b*IT_0027*IT_0070*conj(U_sb_01) +
       mu_h*conj(U_sb_11))));
    const ccomplex_t IT_0210 = -IT_0209;
    const ccomplex_t IT_0211 = IT_0207*IT_0210;
    const ccomplex_t IT_0212 = IT_0205*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u1
      *e_em*IT_0003*IT_0029*IT_0071;
    const ccomplex_t IT_0215 = 0.5*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0003*IT_0029*(mu_h*conj(U_sb_01) + m_b*IT_0027*IT_0071*conj(U_sb_11));
    const ccomplex_t IT_0217 = 0.5*IT_0216;
    const ccomplex_t IT_0218 = IT_0215*IT_0217;
    const ccomplex_t IT_0219 = cpow((-2)*s_13 + IT_0036 + IT_0068 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0220 = IT_0218*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = conj(N_B1)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0223 = IT_0009*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = m_b*conj(N_d1)*e_em*IT_0029*conj(U_sb_01);
    const ccomplex_t IT_0226 = IT_0028*IT_0225;
    const ccomplex_t IT_0227 = 1.4142135623731*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*(IT_0224 + 1.5*IT_0227);
    const ccomplex_t IT_0229 = (-0.333333333333333)*IT_0228;
    const ccomplex_t IT_0230 = IT_0020*IT_0229;
    const ccomplex_t IT_0231 = IT_0038*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = m_b*IT_0232;
    const ccomplex_t IT_0234 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0235 = IT_0009*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = m_b*conj(N_d3)*e_em*IT_0029*conj(U_sb_01);
    const ccomplex_t IT_0238 = IT_0028*IT_0237;
    const ccomplex_t IT_0239 = 1.4142135623731*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*(IT_0236 + 1.5*IT_0239);
    const ccomplex_t IT_0241 = (-0.333333333333333)*IT_0240;
    const ccomplex_t IT_0242 = IT_0053*IT_0241;
    const ccomplex_t IT_0243 = IT_0066*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = m_b*IT_0244;
    const ccomplex_t IT_0246 = IT_0039 + IT_0067 + -IT_0088 + -IT_0105 + 
      -IT_0131 + -IT_0157 + IT_0167 + IT_0181 + IT_0195 + IT_0204 + IT_0213 +
       IT_0221 + -IT_0233 + -IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0003*IT_0027*IT_0029*conj(U_st_00);
    const ccomplex_t IT_0248 = (-0.5)*IT_0247;
    const ccomplex_t IT_0249 = IT_0073*IT_0248;
    const ccomplex_t IT_0250 = IT_0069*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = IT_0009*IT_0070;
    const ccomplex_t IT_0253 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0254 = IT_0252*IT_0253;
    const ccomplex_t IT_0255 = 1.4142135623731*IT_0254;
    const ccomplex_t IT_0256 = IT_0003*IT_0070;
    const ccomplex_t IT_0257 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0258 = IT_0256*IT_0257;
    const ccomplex_t IT_0259 = 1.4142135623731*IT_0258;
    const ccomplex_t IT_0260 = IT_0255 + IT_0259;
    const ccomplex_t IT_0261 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0262 = IT_0256*IT_0261;
    const ccomplex_t IT_0263 = (-2)*IT_0262;
    const ccomplex_t IT_0264 = IT_0260 + IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = 0.5*IT_0265;
    const ccomplex_t IT_0267 = IT_0034*IT_0266;
    const ccomplex_t IT_0268 = IT_0038*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*IT_0268;
    const ccomplex_t IT_0270 = m_b*IT_0269;
    const ccomplex_t IT_0271 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0272 = IT_0256*IT_0271;
    const ccomplex_t IT_0273 = (-2)*IT_0272;
    const ccomplex_t IT_0274 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0275 = IT_0252*IT_0274;
    const ccomplex_t IT_0276 = 1.4142135623731*IT_0275;
    const ccomplex_t IT_0277 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0278 = IT_0256*IT_0277;
    const ccomplex_t IT_0279 = 1.4142135623731*IT_0278;
    const ccomplex_t IT_0280 = IT_0276 + IT_0279;
    const ccomplex_t IT_0281 = IT_0273 + IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = 0.5*IT_0282;
    const ccomplex_t IT_0284 = IT_0179*IT_0283;
    const ccomplex_t IT_0285 = IT_0106*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = m_b*IT_0286;
    const ccomplex_t IT_0288 = IT_0229*IT_0266;
    const ccomplex_t IT_0289 = IT_0000*IT_0038*IT_0288;
    const ccomplex_t IT_0290 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0291 = IT_0256*IT_0290;
    const ccomplex_t IT_0292 = (-2)*IT_0291;
    const ccomplex_t IT_0293 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0294 = IT_0252*IT_0293;
    const ccomplex_t IT_0295 = 1.4142135623731*IT_0294;
    const ccomplex_t IT_0296 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0297 = IT_0256*IT_0296;
    const ccomplex_t IT_0298 = 1.4142135623731*IT_0297;
    const ccomplex_t IT_0299 = IT_0295 + IT_0298;
    const ccomplex_t IT_0300 = IT_0292 + IT_0299;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*IT_0300;
    const ccomplex_t IT_0302 = 0.5*IT_0301;
    const ccomplex_t IT_0303 = IT_0152*IT_0302;
    const ccomplex_t IT_0304 = IT_0154*IT_0182*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_cb*e_em*IT_0003*IT_0027*IT_0029;
    const ccomplex_t IT_0306 = (-0.5)*IT_0305;
    const ccomplex_t IT_0307 = IT_0210*IT_0306;
    const ccomplex_t IT_0308 = IT_0205*IT_0307;
    const ccomplex_t IT_0309 = (0 + _Complex_I*1)*IT_0308;
    const ccomplex_t IT_0310 = IT_0114*IT_0283;
    const ccomplex_t IT_0311 = IT_0106*IT_0168*IT_0310;
    const ccomplex_t IT_0312 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0313 = IT_0256*IT_0312;
    const ccomplex_t IT_0314 = 1.4142135623731*IT_0313;
    const ccomplex_t IT_0315 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0316 = IT_0252*IT_0315;
    const ccomplex_t IT_0317 = 1.4142135623731*IT_0316;
    const ccomplex_t IT_0318 = IT_0314 + IT_0317;
    const ccomplex_t IT_0319 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0320 = IT_0256*IT_0319;
    const ccomplex_t IT_0321 = (-2)*IT_0320;
    const ccomplex_t IT_0322 = IT_0318 + IT_0321;
    const ccomplex_t IT_0323 = (0 + _Complex_I*1)*IT_0322;
    const ccomplex_t IT_0324 = 0.5*IT_0323;
    const ccomplex_t IT_0325 = IT_0241*IT_0324;
    const ccomplex_t IT_0326 = IT_0040*IT_0066*IT_0325;
    const ccomplex_t IT_0327 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *e_em*IT_0003*IT_0027*IT_0029*IT_0158*V_ub_mod;
    const ccomplex_t IT_0328 = (-0.5)*IT_0327;
    const ccomplex_t IT_0329 = IT_0201*IT_0328;
    const ccomplex_t IT_0330 = IT_0196*IT_0329;
    const ccomplex_t IT_0331 = (0 + _Complex_I*1)*IT_0330;
    const ccomplex_t IT_0332 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0003*IT_0027*IT_0029*conj(U_st_01);
    const ccomplex_t IT_0333 = (-0.5)*IT_0332;
    const ccomplex_t IT_0334 = IT_0091*IT_0333;
    const ccomplex_t IT_0335 = IT_0089*IT_0334;
    const ccomplex_t IT_0336 = (0 + _Complex_I*1)*IT_0335;
    const ccomplex_t IT_0337 = IT_0064*IT_0324;
    const ccomplex_t IT_0338 = IT_0066*IT_0337;
    const ccomplex_t IT_0339 = (0 + _Complex_I*1)*IT_0338;
    const ccomplex_t IT_0340 = m_b*IT_0339;
    const ccomplex_t IT_0341 = IT_0193*IT_0302;
    const ccomplex_t IT_0342 = IT_0154*IT_0341;
    const ccomplex_t IT_0343 = (0 + _Complex_I*1)*IT_0342;
    const ccomplex_t IT_0344 = m_b*IT_0343;
    const ccomplex_t IT_0345 = -IT_0251 + -IT_0270 + -IT_0287 + IT_0289 +
       IT_0304 + -IT_0309 + IT_0311 + IT_0326 + -IT_0331 + -IT_0336 + -IT_0340 +
       -IT_0344;
    const ccomplex_t IT_0346 = m_b*m_C_1;
    const ccomplex_t IT_0347 = -IT_0269 + -IT_0286 + -IT_0339 + -IT_0343;
    const ccomplex_t IT_0348 = s_23*m_C_1;
    const ccomplex_t IT_0349 = -IT_0130 + -IT_0156 + -IT_0232 + -IT_0244;
    const ccomplex_t IT_0350 = m_b*s_12;
    const ccomplex_t IT_0351 = 6*s_13;
    const ccomplex_t IT_0352 = (-6)*IT_0345;
    const ccomplex_t IT_0353 = m_b*m_C_1*IT_0037;
    const ccomplex_t IT_0354 = s_13*IT_0037;
    const ccomplex_t IT_0355 = (-6)*IT_0354;
    const ccomplex_t IT_0356 = s_12*s_23;
    const ccomplex_t IT_0357 = 12*IT_0356;
    const ccomplex_t IT_0358 = IT_0355 + IT_0357;
    const ccomplex_t IT_0359 = (-6)*conj(IT_0345);
    const ccomplex_t IT_0360 = 0.166666666666667*IT_0348;
    const ccomplex_t IT_0361 = (-6)*conj(IT_0246);
    const ccomplex_t IT_0362 = 0.166666666666667*IT_0350;
    const ccomplex_t IT_0363 = 6*IT_0246*(conj(IT_0345)*IT_0346 + -conj
      (IT_0347)*IT_0348 + -conj(IT_0349)*IT_0350 + 0.166666666666667*conj
      (IT_0246)*IT_0351) + 6*IT_0345*(conj(IT_0246)*IT_0346 + 0.166666666666667
      *conj(IT_0345)*IT_0351) + (IT_0348*conj(IT_0349) + conj(IT_0347)*IT_0350)
      *IT_0352 + 6*IT_0347*(conj(IT_0349)*IT_0353 + 0.166666666666667*conj
      (IT_0347)*IT_0358 + IT_0360*IT_0361 + IT_0359*IT_0362) + 6*IT_0349*(conj
      (IT_0347)*IT_0353 + 0.166666666666667*conj(IT_0349)*IT_0358 + IT_0359
      *IT_0360 + IT_0361*IT_0362);
    return create_ccomplex_return(IT_0363);
}

