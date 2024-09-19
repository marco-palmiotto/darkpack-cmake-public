#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_ss_R_to_anti_s_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_ss_R_to_anti_s_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
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
    const creal_t m_ss_R = param->m_ss_R;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp2*IT_0002;
    const ccomplex_t IT_0006 = -IT_0005;
    const ccomplex_t IT_0007 = IT_0004*IT_0006;
    const ccomplex_t IT_0008 = pow(m_s, 2);
    const ccomplex_t IT_0009 = pow(m_C_2, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_13 + IT_0008 + IT_0009 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = sin(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_c*m_s*conj
      (V_cs)*e_em*IT_0000*IT_0002*IT_0014*IT_0016;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0000*IT_0002*IT_0016;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0008 + IT_0009 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_s*m_u*V_us
      *e_em*IT_0000*IT_0002*IT_0014*IT_0016;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_us
      *e_em*IT_0000*IT_0002*IT_0016;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0008 + IT_0009 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0000*IT_0002*IT_0014;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0002*IT_0013;
    const ccomplex_t IT_0037 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = 2*IT_0038;
    const ccomplex_t IT_0040 = cos(theta_W);
    const ccomplex_t IT_0041 = cpow(IT_0040, -1);
    const ccomplex_t IT_0042 = IT_0013*IT_0041;
    const ccomplex_t IT_0043 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0047 = IT_0036*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = IT_0045 + IT_0048;
    const ccomplex_t IT_0050 = IT_0039 + IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = IT_0035*IT_0052;
    const ccomplex_t IT_0054 = pow(m_ss_R, 2);
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0008 + IT_0054 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0033*IT_0053*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0008 + IT_0054 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0000*IT_0002*IT_0014;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0062 = IT_0036*IT_0061;
    const ccomplex_t IT_0063 = 2*IT_0062;
    const ccomplex_t IT_0064 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0065 = IT_0042*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0068 = IT_0036*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = IT_0066 + IT_0069;
    const ccomplex_t IT_0071 = IT_0063 + IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = IT_0060*IT_0073;
    const ccomplex_t IT_0075 = IT_0057*IT_0058*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0041;
    const ccomplex_t IT_0077 = (-0.333333333333333)*IT_0076;
    const ccomplex_t IT_0078 = IT_0073*IT_0077;
    const ccomplex_t IT_0079 = IT_0058*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = m_s*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_13 + IT_0008 + IT_0009 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *mu_h*IT_0000*IT_0002;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*V_us*e_em*V_Wp2*IT_0002;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = IT_0082*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0091 = cpow((-2)*s_23 + IT_0008 + IT_0054 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0000*IT_0002*IT_0014;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0095 = IT_0036*IT_0094;
    const ccomplex_t IT_0096 = 2*IT_0095;
    const ccomplex_t IT_0097 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0098 = IT_0042*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0101 = IT_0036*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = IT_0099 + IT_0102;
    const ccomplex_t IT_0104 = IT_0096 + IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0093*IT_0106;
    const ccomplex_t IT_0108 = IT_0090*IT_0091*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0110 = cpow((-2)*s_23 + IT_0008 + IT_0054 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0000*IT_0002*IT_0014;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0114 = IT_0036*IT_0113;
    const ccomplex_t IT_0115 = 2*IT_0114;
    const ccomplex_t IT_0116 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0117 = IT_0042*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0120 = IT_0036*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = IT_0118 + IT_0121;
    const ccomplex_t IT_0123 = IT_0115 + IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = IT_0112*IT_0125;
    const ccomplex_t IT_0127 = IT_0109*IT_0110*IT_0126;
    const ccomplex_t IT_0128 = cpow((-2)*s_13 + IT_0008 + IT_0009 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*conj(V_ts)*(A_b*IT_0015
      *U_st_00 + 0.707106781186548*m_s*e_em*IT_0000*IT_0002*(mu_h*U_st_00 + m_t
      *IT_0014*IT_0016*U_st_10));
    const ccomplex_t IT_0130 = IT_0002*IT_0016;
    const ccomplex_t IT_0131 = m_t*V_ts*V_u2*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0132 = IT_0130*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = V_ts*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0137 = IT_0002*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = -IT_0138;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = IT_0135 + IT_0140;
    const ccomplex_t IT_0142 = IT_0129*IT_0141;
    const ccomplex_t IT_0143 = IT_0128*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = cpow((-2)*s_13 + IT_0008 + IT_0009 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*conj(V_ts)*(A_b*IT_0015
      *U_st_01 + 0.707106781186548*m_s*e_em*IT_0000*IT_0002*(mu_h*U_st_01 + m_t
      *IT_0014*IT_0016*U_st_11));
    const ccomplex_t IT_0147 = m_t*V_ts*V_u2*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0148 = IT_0130*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = V_ts*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0153 = IT_0002*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0151 + IT_0156;
    const ccomplex_t IT_0158 = IT_0146*IT_0157;
    const ccomplex_t IT_0159 = IT_0145*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0041;
    const ccomplex_t IT_0162 = (-0.333333333333333)*IT_0161;
    const ccomplex_t IT_0163 = IT_0052*IT_0162;
    const ccomplex_t IT_0164 = IT_0055*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = m_s*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0041;
    const ccomplex_t IT_0168 = (-0.333333333333333)*IT_0167;
    const ccomplex_t IT_0169 = IT_0106*IT_0168;
    const ccomplex_t IT_0170 = IT_0091*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = m_s*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0041;
    const ccomplex_t IT_0174 = (-0.333333333333333)*IT_0173;
    const ccomplex_t IT_0175 = IT_0125*IT_0174;
    const ccomplex_t IT_0176 = IT_0110*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = m_s*IT_0177;
    const ccomplex_t IT_0179 = IT_0012 + IT_0024 + IT_0032 + IT_0056 + IT_0075
       + -IT_0081 + IT_0089 + IT_0108 + IT_0127 + -IT_0144 + -IT_0160 + -IT_0166
       + -IT_0172 + -IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0000*IT_0002*IT_0014;
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = IT_0004*IT_0181;
    const ccomplex_t IT_0183 = IT_0010*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_us*e_em*IT_0000*IT_0002*IT_0014;
    const ccomplex_t IT_0186 = (-0.5)*IT_0185;
    const ccomplex_t IT_0187 = IT_0084*IT_0186;
    const ccomplex_t IT_0188 = IT_0082*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = IT_0002*IT_0015;
    const ccomplex_t IT_0191 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0192 = IT_0190*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = IT_0015*IT_0041;
    const ccomplex_t IT_0195 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0196 = IT_0194*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = IT_0193 + IT_0197;
    const ccomplex_t IT_0199 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0200 = IT_0190*IT_0199;
    const ccomplex_t IT_0201 = (-2)*IT_0200;
    const ccomplex_t IT_0202 = IT_0198 + IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = 0.5*IT_0203;
    const ccomplex_t IT_0205 = IT_0162*IT_0204;
    const ccomplex_t IT_0206 = IT_0033*IT_0055*IT_0205;
    const ccomplex_t IT_0207 = IT_0035*IT_0204;
    const ccomplex_t IT_0208 = IT_0055*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = m_s*IT_0209;
    const ccomplex_t IT_0211 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0212 = IT_0190*IT_0211;
    const ccomplex_t IT_0213 = (-2)*IT_0212;
    const ccomplex_t IT_0214 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0215 = IT_0194*IT_0214;
    const ccomplex_t IT_0216 = 1.4142135623731*IT_0215;
    const ccomplex_t IT_0217 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0218 = IT_0190*IT_0217;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = IT_0216 + IT_0219;
    const ccomplex_t IT_0221 = IT_0213 + IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = 0.5*IT_0222;
    const ccomplex_t IT_0224 = IT_0174*IT_0223;
    const ccomplex_t IT_0225 = IT_0109*IT_0110*IT_0224;
    const ccomplex_t IT_0226 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0227 = IT_0190*IT_0226;
    const ccomplex_t IT_0228 = (-2)*IT_0227;
    const ccomplex_t IT_0229 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0230 = IT_0194*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0233 = IT_0190*IT_0232;
    const ccomplex_t IT_0234 = 1.4142135623731*IT_0233;
    const ccomplex_t IT_0235 = IT_0231 + IT_0234;
    const ccomplex_t IT_0236 = IT_0228 + IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = 0.5*IT_0237;
    const ccomplex_t IT_0239 = IT_0077*IT_0238;
    const ccomplex_t IT_0240 = IT_0057*IT_0058*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0000*IT_0002*IT_0014*conj(U_st_00);
    const ccomplex_t IT_0242 = (-0.5)*IT_0241;
    const ccomplex_t IT_0243 = IT_0129*IT_0242;
    const ccomplex_t IT_0244 = IT_0128*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0000*IT_0002*IT_0014*conj(U_st_01);
    const ccomplex_t IT_0247 = (-0.5)*IT_0246;
    const ccomplex_t IT_0248 = IT_0146*IT_0247;
    const ccomplex_t IT_0249 = IT_0145*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*IT_0249;
    const ccomplex_t IT_0251 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0252 = IT_0190*IT_0251;
    const ccomplex_t IT_0253 = (-2)*IT_0252;
    const ccomplex_t IT_0254 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0255 = IT_0194*IT_0254;
    const ccomplex_t IT_0256 = 1.4142135623731*IT_0255;
    const ccomplex_t IT_0257 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0258 = IT_0190*IT_0257;
    const ccomplex_t IT_0259 = 1.4142135623731*IT_0258;
    const ccomplex_t IT_0260 = IT_0256 + IT_0259;
    const ccomplex_t IT_0261 = IT_0253 + IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = 0.5*IT_0262;
    const ccomplex_t IT_0264 = IT_0093*IT_0263;
    const ccomplex_t IT_0265 = IT_0091*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = m_s*IT_0266;
    const ccomplex_t IT_0268 = IT_0112*IT_0223;
    const ccomplex_t IT_0269 = IT_0110*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*IT_0269;
    const ccomplex_t IT_0271 = m_s*IT_0270;
    const ccomplex_t IT_0272 = IT_0060*IT_0238;
    const ccomplex_t IT_0273 = IT_0058*IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*IT_0273;
    const ccomplex_t IT_0275 = m_s*IT_0274;
    const ccomplex_t IT_0276 = IT_0168*IT_0263;
    const ccomplex_t IT_0277 = IT_0090*IT_0091*IT_0276;
    const ccomplex_t IT_0278 = -IT_0184 + -IT_0189 + IT_0206 + -IT_0210 +
       IT_0225 + IT_0240 + -IT_0245 + -IT_0250 + -IT_0267 + -IT_0271 + -IT_0275 
      + IT_0277;
    const ccomplex_t IT_0279 = m_s*m_C_2;
    const ccomplex_t IT_0280 = -IT_0209 + -IT_0266 + -IT_0270 + -IT_0274;
    const ccomplex_t IT_0281 = s_23*m_C_2;
    const ccomplex_t IT_0282 = -IT_0080 + -IT_0165 + -IT_0171 + -IT_0177;
    const ccomplex_t IT_0283 = m_s*s_12;
    const ccomplex_t IT_0284 = 6*s_13;
    const ccomplex_t IT_0285 = (-6)*IT_0278;
    const ccomplex_t IT_0286 = m_s*m_C_2*IT_0054;
    const ccomplex_t IT_0287 = s_13*IT_0054;
    const ccomplex_t IT_0288 = (-6)*IT_0287;
    const ccomplex_t IT_0289 = s_12*s_23;
    const ccomplex_t IT_0290 = 12*IT_0289;
    const ccomplex_t IT_0291 = IT_0288 + IT_0290;
    const ccomplex_t IT_0292 = (-6)*conj(IT_0278);
    const ccomplex_t IT_0293 = 0.166666666666667*IT_0281;
    const ccomplex_t IT_0294 = (-6)*conj(IT_0179);
    const ccomplex_t IT_0295 = 0.166666666666667*IT_0283;
    const ccomplex_t IT_0296 = 6*IT_0179*(conj(IT_0278)*IT_0279 + -conj
      (IT_0280)*IT_0281 + -conj(IT_0282)*IT_0283 + 0.166666666666667*conj
      (IT_0179)*IT_0284) + 6*IT_0278*(conj(IT_0179)*IT_0279 + 0.166666666666667
      *conj(IT_0278)*IT_0284) + (IT_0281*conj(IT_0282) + conj(IT_0280)*IT_0283)
      *IT_0285 + 6*IT_0280*(conj(IT_0282)*IT_0286 + 0.166666666666667*conj
      (IT_0280)*IT_0291 + IT_0293*IT_0294 + IT_0292*IT_0295) + 6*IT_0282*(conj
      (IT_0280)*IT_0286 + 0.166666666666667*conj(IT_0282)*IT_0291 + IT_0292
      *IT_0293 + IT_0294*IT_0295);
    return create_ccomplex_return(IT_0296);
}

