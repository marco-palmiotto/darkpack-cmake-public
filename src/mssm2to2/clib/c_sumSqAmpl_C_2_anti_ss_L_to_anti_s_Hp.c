#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_ss_L_to_anti_s_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_ss_L_to_anti_s_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_L = param->m_ss_L;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0001 = pow(m_s, 2);
    const ccomplex_t IT_0002 = pow(m_ss_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = N_B4*e_em;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = N_W4*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0008 + (-3)*IT_0013);
    const ccomplex_t IT_0015 = 0.166666666666667*IT_0014;
    const ccomplex_t IT_0016 = cos(beta);
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 2*IT_0019;
    const ccomplex_t IT_0021 = IT_0005*IT_0016;
    const ccomplex_t IT_0022 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0026 = IT_0017*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = IT_0024 + IT_0027;
    const ccomplex_t IT_0029 = IT_0020 + IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = IT_0015*IT_0031;
    const ccomplex_t IT_0033 = IT_0000*IT_0003*IT_0032;
    const ccomplex_t IT_0034 = pow(m_W, -1);
    const ccomplex_t IT_0035 = cpow(IT_0016, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0010*IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0031*IT_0037;
    const ccomplex_t IT_0039 = IT_0003*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = m_s*IT_0040;
    const ccomplex_t IT_0042 = sin(beta);
    const ccomplex_t IT_0043 = pow(m_c, 2);
    const ccomplex_t IT_0044 = cpow(IT_0042, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0010*(IT_0001*IT_0034*IT_0035*IT_0042 + (-2)*IT_0016*(m_W*IT_0042
       + (-0.5)*IT_0034*IT_0043*IT_0044));
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp2*IT_0010;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = pow(m_C_2, 2);
    const ccomplex_t IT_0051 = cpow((-2)*s_13 + IT_0001 + IT_0050 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = pow(m_u, 2);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0010*(IT_0001*IT_0034*IT_0035*IT_0042 + (-2)*IT_0016*(m_W*IT_0042 + (
      -0.5)*IT_0034*IT_0044*IT_0054));
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*V_us*e_em*V_Wp2*IT_0010;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = IT_0056*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0001 + IT_0050 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*mu_h*IT_0010*IT_0034;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0010*IT_0034*IT_0044;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0064*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_13 + IT_0001 + IT_0050 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *mu_h*IT_0010*IT_0034;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_us
      *e_em*IT_0010*IT_0034*IT_0044;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_13 + IT_0001 + IT_0050 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0080 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0081 = N_B1*e_em;
    const ccomplex_t IT_0082 = IT_0005*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = N_W1*e_em;
    const ccomplex_t IT_0085 = IT_0010*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + (-3)*IT_0086);
    const ccomplex_t IT_0088 = 0.166666666666667*IT_0087;
    const ccomplex_t IT_0089 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0090 = IT_0017*IT_0089;
    const ccomplex_t IT_0091 = 2*IT_0090;
    const ccomplex_t IT_0092 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0093 = IT_0021*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0096 = IT_0017*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = IT_0094 + IT_0097;
    const ccomplex_t IT_0099 = IT_0091 + IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = IT_0088*IT_0101;
    const ccomplex_t IT_0103 = IT_0079*IT_0080*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0105 = N_B3*e_em;
    const ccomplex_t IT_0106 = IT_0005*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = N_W3*e_em;
    const ccomplex_t IT_0109 = IT_0010*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0107 + (-3)*IT_0110);
    const ccomplex_t IT_0112 = 0.166666666666667*IT_0111;
    const ccomplex_t IT_0113 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0114 = IT_0017*IT_0113;
    const ccomplex_t IT_0115 = 2*IT_0114;
    const ccomplex_t IT_0116 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0117 = IT_0021*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0120 = IT_0017*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = IT_0118 + IT_0121;
    const ccomplex_t IT_0123 = IT_0115 + IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = IT_0112*IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0128 = IT_0104*IT_0126*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0130 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0131 = N_B2*e_em;
    const ccomplex_t IT_0132 = IT_0005*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = N_W2*e_em;
    const ccomplex_t IT_0135 = IT_0010*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + (-3)*IT_0136);
    const ccomplex_t IT_0138 = 0.166666666666667*IT_0137;
    const ccomplex_t IT_0139 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0140 = IT_0017*IT_0139;
    const ccomplex_t IT_0141 = 2*IT_0140;
    const ccomplex_t IT_0142 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0143 = IT_0021*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0146 = IT_0017*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = IT_0144 + IT_0147;
    const ccomplex_t IT_0149 = IT_0141 + IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0138*IT_0151;
    const ccomplex_t IT_0153 = IT_0129*IT_0130*IT_0152;
    const ccomplex_t IT_0154 = cpow((-2)*s_13 + IT_0001 + IT_0050 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0155 = A_t*IT_0016*U_st_10 + (-1.4142135623731)*e_em
      *IT_0010*(m_W*IT_0016*IT_0042*U_st_00 + (-0.5)*IT_0034*(IT_0001*IT_0035
      *IT_0042*U_st_00 + m_t*(m_t*IT_0016*IT_0044*U_st_00 + mu_h*U_st_10)));
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*conj(V_ts)*IT_0155;
    const ccomplex_t IT_0157 = IT_0010*IT_0044;
    const ccomplex_t IT_0158 = m_t*V_ts*V_u2*e_em*IT_0034*conj(U_st_10);
    const ccomplex_t IT_0159 = IT_0157*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = V_ts*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0164 = IT_0010*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0162 + IT_0167;
    const ccomplex_t IT_0169 = IT_0156*IT_0168;
    const ccomplex_t IT_0170 = IT_0154*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = cpow((-2)*s_13 + IT_0001 + IT_0050 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0173 = A_t*IT_0016*U_st_11 + (-1.4142135623731)*e_em
      *IT_0010*(m_W*IT_0016*IT_0042*U_st_01 + (-0.5)*IT_0034*(IT_0001*IT_0035
      *IT_0042*U_st_01 + m_t*(m_t*IT_0016*IT_0044*U_st_01 + mu_h*U_st_11)));
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*conj(V_ts)*IT_0173;
    const ccomplex_t IT_0175 = m_t*V_ts*V_u2*e_em*IT_0034*conj(U_st_11);
    const ccomplex_t IT_0176 = IT_0157*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = (-0.5)*IT_0178;
    const ccomplex_t IT_0180 = V_ts*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0181 = IT_0010*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = -IT_0182;
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = IT_0179 + IT_0184;
    const ccomplex_t IT_0186 = IT_0174*IT_0185;
    const ccomplex_t IT_0187 = IT_0172*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0010*IT_0034*IT_0035;
    const ccomplex_t IT_0190 = (-0.5)*IT_0189;
    const ccomplex_t IT_0191 = IT_0101*IT_0190;
    const ccomplex_t IT_0192 = IT_0080*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = m_s*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0010*IT_0034*IT_0035;
    const ccomplex_t IT_0196 = (-0.5)*IT_0195;
    const ccomplex_t IT_0197 = IT_0151*IT_0196;
    const ccomplex_t IT_0198 = IT_0130*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = m_s*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0010*IT_0034*IT_0035;
    const ccomplex_t IT_0202 = (-0.5)*IT_0201;
    const ccomplex_t IT_0203 = IT_0125*IT_0202;
    const ccomplex_t IT_0204 = IT_0127*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = m_s*IT_0205;
    const ccomplex_t IT_0207 = IT_0033 + -IT_0041 + IT_0053 + IT_0062 +
       IT_0070 + IT_0078 + IT_0103 + IT_0128 + IT_0153 + -IT_0171 + -IT_0188 + 
      -IT_0194 + -IT_0200 + -IT_0206;
    const ccomplex_t IT_0208 = IT_0010*IT_0042;
    const ccomplex_t IT_0209 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0210 = IT_0208*IT_0209;
    const ccomplex_t IT_0211 = (-2)*IT_0210;
    const ccomplex_t IT_0212 = IT_0005*IT_0042;
    const ccomplex_t IT_0213 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0214 = IT_0212*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0217 = IT_0208*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = IT_0215 + IT_0218;
    const ccomplex_t IT_0220 = IT_0211 + IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = 0.5*IT_0221;
    const ccomplex_t IT_0223 = IT_0037*IT_0222;
    const ccomplex_t IT_0224 = IT_0000*IT_0003*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0010*IT_0034*IT_0035*conj(U_st_00);
    const ccomplex_t IT_0226 = (-0.5)*IT_0225;
    const ccomplex_t IT_0227 = IT_0156*IT_0226;
    const ccomplex_t IT_0228 = IT_0154*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0010*IT_0034*IT_0035*conj(U_st_01);
    const ccomplex_t IT_0231 = (-0.5)*IT_0230;
    const ccomplex_t IT_0232 = IT_0174*IT_0231;
    const ccomplex_t IT_0233 = IT_0172*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0010*IT_0034*IT_0035;
    const ccomplex_t IT_0236 = (-0.5)*IT_0235;
    const ccomplex_t IT_0237 = IT_0046*IT_0236;
    const ccomplex_t IT_0238 = IT_0051*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_us*e_em*IT_0010*IT_0034*IT_0035;
    const ccomplex_t IT_0241 = (-0.5)*IT_0240;
    const ccomplex_t IT_0242 = IT_0056*IT_0241;
    const ccomplex_t IT_0243 = IT_0060*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0246 = IT_0208*IT_0245;
    const ccomplex_t IT_0247 = 1.4142135623731*IT_0246;
    const ccomplex_t IT_0248 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0249 = IT_0212*IT_0248;
    const ccomplex_t IT_0250 = 1.4142135623731*IT_0249;
    const ccomplex_t IT_0251 = IT_0247 + IT_0250;
    const ccomplex_t IT_0252 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0253 = IT_0208*IT_0252;
    const ccomplex_t IT_0254 = (-2)*IT_0253;
    const ccomplex_t IT_0255 = IT_0251 + IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = 0.5*IT_0256;
    const ccomplex_t IT_0258 = IT_0190*IT_0257;
    const ccomplex_t IT_0259 = IT_0079*IT_0080*IT_0258;
    const ccomplex_t IT_0260 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0261 = IT_0208*IT_0260;
    const ccomplex_t IT_0262 = (-2)*IT_0261;
    const ccomplex_t IT_0263 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0264 = IT_0212*IT_0263;
    const ccomplex_t IT_0265 = 1.4142135623731*IT_0264;
    const ccomplex_t IT_0266 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0267 = IT_0208*IT_0266;
    const ccomplex_t IT_0268 = 1.4142135623731*IT_0267;
    const ccomplex_t IT_0269 = IT_0265 + IT_0268;
    const ccomplex_t IT_0270 = IT_0262 + IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = 0.5*IT_0271;
    const ccomplex_t IT_0273 = IT_0196*IT_0272;
    const ccomplex_t IT_0274 = IT_0129*IT_0130*IT_0273;
    const ccomplex_t IT_0275 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0276 = IT_0208*IT_0275;
    const ccomplex_t IT_0277 = (-2)*IT_0276;
    const ccomplex_t IT_0278 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0279 = IT_0212*IT_0278;
    const ccomplex_t IT_0280 = 1.4142135623731*IT_0279;
    const ccomplex_t IT_0281 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0282 = IT_0208*IT_0281;
    const ccomplex_t IT_0283 = 1.4142135623731*IT_0282;
    const ccomplex_t IT_0284 = IT_0280 + IT_0283;
    const ccomplex_t IT_0285 = IT_0277 + IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = 0.5*IT_0286;
    const ccomplex_t IT_0288 = IT_0202*IT_0287;
    const ccomplex_t IT_0289 = IT_0104*IT_0127*IT_0288;
    const ccomplex_t IT_0290 = IT_0088*IT_0257;
    const ccomplex_t IT_0291 = IT_0080*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = m_s*IT_0292;
    const ccomplex_t IT_0294 = IT_0138*IT_0272;
    const ccomplex_t IT_0295 = IT_0130*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = m_s*IT_0296;
    const ccomplex_t IT_0298 = IT_0112*IT_0287;
    const ccomplex_t IT_0299 = IT_0127*IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = m_s*IT_0300;
    const ccomplex_t IT_0302 = IT_0015*IT_0222;
    const ccomplex_t IT_0303 = IT_0003*IT_0302;
    const ccomplex_t IT_0304 = (0 + _Complex_I*1)*IT_0303;
    const ccomplex_t IT_0305 = m_s*IT_0304;
    const ccomplex_t IT_0306 = IT_0224 + -IT_0229 + -IT_0234 + -IT_0239 + 
      -IT_0244 + IT_0259 + IT_0274 + IT_0289 + -IT_0293 + -IT_0297 + -IT_0301 + 
      -IT_0305;
    const ccomplex_t IT_0307 = m_s*m_C_2;
    const ccomplex_t IT_0308 = -IT_0292 + -IT_0296 + -IT_0300 + -IT_0304;
    const ccomplex_t IT_0309 = s_23*m_C_2;
    const ccomplex_t IT_0310 = -IT_0040 + -IT_0193 + -IT_0199 + -IT_0205;
    const ccomplex_t IT_0311 = m_s*s_12;
    const ccomplex_t IT_0312 = 6*s_13;
    const ccomplex_t IT_0313 = (-6)*IT_0306;
    const ccomplex_t IT_0314 = m_s*m_C_2*IT_0002;
    const ccomplex_t IT_0315 = s_13*IT_0002;
    const ccomplex_t IT_0316 = (-6)*IT_0315;
    const ccomplex_t IT_0317 = s_12*s_23;
    const ccomplex_t IT_0318 = 12*IT_0317;
    const ccomplex_t IT_0319 = IT_0316 + IT_0318;
    const ccomplex_t IT_0320 = (-6)*conj(IT_0306);
    const ccomplex_t IT_0321 = 0.166666666666667*IT_0309;
    const ccomplex_t IT_0322 = (-6)*conj(IT_0207);
    const ccomplex_t IT_0323 = 0.166666666666667*IT_0311;
    const ccomplex_t IT_0324 = 6*IT_0207*(conj(IT_0306)*IT_0307 + -conj
      (IT_0308)*IT_0309 + -conj(IT_0310)*IT_0311 + 0.166666666666667*conj
      (IT_0207)*IT_0312) + 6*IT_0306*(conj(IT_0207)*IT_0307 + 0.166666666666667
      *conj(IT_0306)*IT_0312) + (IT_0309*conj(IT_0310) + conj(IT_0308)*IT_0311)
      *IT_0313 + 6*IT_0308*(conj(IT_0310)*IT_0314 + 0.166666666666667*conj
      (IT_0308)*IT_0319 + IT_0321*IT_0322 + IT_0320*IT_0323) + 6*IT_0310*(conj
      (IT_0308)*IT_0314 + 0.166666666666667*conj(IT_0310)*IT_0319 + IT_0320
      *IT_0321 + IT_0322*IT_0323);
    return create_ccomplex_return(IT_0324);
}

