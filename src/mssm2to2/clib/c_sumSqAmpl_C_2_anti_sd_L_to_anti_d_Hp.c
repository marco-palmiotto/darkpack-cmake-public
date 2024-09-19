#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_sd_L_to_anti_d_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_sd_L_to_anti_d_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
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
    const creal_t m_sd_L = param->m_sd_L;
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
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = pow(m_d, 2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = pow(m_c, 2);
    const ccomplex_t IT_0006 = cpow(IT_0004, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*conj(V_cd)
      *e_em*(IT_0000*IT_0001*IT_0003*IT_0004 + (-2)*IT_0002*(m_W*IT_0004 + (-0.5
      )*IT_0000*IT_0005*IT_0006))*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*V_cd*e_em*V_Wp2*IT_0008;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = pow(m_C_2, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0001 + IT_0014 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *e_em*mu_h*IT_0000*IT_0008;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u2
      *e_em*IT_0000*IT_0006*IT_0008;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0001 + IT_0014 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *mu_h*IT_0000*IT_0008;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0000*IT_0006*IT_0008;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0001 + IT_0014 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = pow(m_u, 2);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0008*(IT_0000*IT_0001*IT_0003*IT_0004 + (-2)*IT_0002*(m_W*IT_0004 + (
      -0.5)*IT_0000*IT_0006*IT_0034));
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0008;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0001 + IT_0014 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0044 = pow(m_sd_L, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0001 + IT_0044 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0046 = cos(theta_W);
    const ccomplex_t IT_0047 = cpow(IT_0046, -1);
    const ccomplex_t IT_0048 = N_B1*e_em;
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = N_W1*e_em;
    const ccomplex_t IT_0052 = IT_0008*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + (-3)*IT_0053);
    const ccomplex_t IT_0055 = 0.166666666666667*IT_0054;
    const ccomplex_t IT_0056 = IT_0002*IT_0008;
    const ccomplex_t IT_0057 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = 2*IT_0058;
    const ccomplex_t IT_0060 = IT_0002*IT_0047;
    const ccomplex_t IT_0061 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0065 = IT_0056*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = IT_0063 + IT_0066;
    const ccomplex_t IT_0068 = IT_0059 + IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = IT_0055*IT_0070;
    const ccomplex_t IT_0072 = IT_0043*IT_0045*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0074 = cpow((-2)*s_23 + IT_0001 + IT_0044 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0075 = N_B3*e_em;
    const ccomplex_t IT_0076 = IT_0047*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = N_W3*e_em;
    const ccomplex_t IT_0079 = IT_0008*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + (-3)*IT_0080);
    const ccomplex_t IT_0082 = 0.166666666666667*IT_0081;
    const ccomplex_t IT_0083 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0084 = IT_0056*IT_0083;
    const ccomplex_t IT_0085 = 2*IT_0084;
    const ccomplex_t IT_0086 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0087 = IT_0060*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0090 = IT_0056*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = IT_0088 + IT_0091;
    const ccomplex_t IT_0093 = IT_0085 + IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = IT_0082*IT_0095;
    const ccomplex_t IT_0097 = IT_0073*IT_0074*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0000*IT_0003*IT_0008;
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0070*IT_0099;
    const ccomplex_t IT_0101 = IT_0045*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = m_d*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0105 = cpow((-2)*s_23 + IT_0001 + IT_0044 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0106 = N_B2*e_em;
    const ccomplex_t IT_0107 = IT_0047*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = N_W2*e_em;
    const ccomplex_t IT_0110 = IT_0008*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + (-3)*IT_0111);
    const ccomplex_t IT_0113 = 0.166666666666667*IT_0112;
    const ccomplex_t IT_0114 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0115 = IT_0056*IT_0114;
    const ccomplex_t IT_0116 = 2*IT_0115;
    const ccomplex_t IT_0117 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0118 = IT_0060*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0121 = IT_0056*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = IT_0119 + IT_0122;
    const ccomplex_t IT_0124 = IT_0116 + IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = 0.5*IT_0125;
    const ccomplex_t IT_0127 = IT_0113*IT_0126;
    const ccomplex_t IT_0128 = IT_0104*IT_0105*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0130 = cpow((-2)*s_23 + IT_0001 + IT_0044 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0131 = N_B4*e_em;
    const ccomplex_t IT_0132 = IT_0047*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = N_W4*e_em;
    const ccomplex_t IT_0135 = IT_0008*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + (-3)*IT_0136);
    const ccomplex_t IT_0138 = 0.166666666666667*IT_0137;
    const ccomplex_t IT_0139 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0140 = IT_0056*IT_0139;
    const ccomplex_t IT_0141 = 2*IT_0140;
    const ccomplex_t IT_0142 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0143 = IT_0060*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0146 = IT_0056*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = IT_0144 + IT_0147;
    const ccomplex_t IT_0149 = IT_0141 + IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0138*IT_0151;
    const ccomplex_t IT_0153 = IT_0129*IT_0130*IT_0152;
    const ccomplex_t IT_0154 = cpow((-2)*s_13 + IT_0001 + IT_0014 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0155 = A_t*IT_0002*U_st_10 + (-1.4142135623731)*e_em
      *IT_0008*(m_W*IT_0002*IT_0004*U_st_00 + (-0.5)*IT_0000*(IT_0001*IT_0003
      *IT_0004*U_st_00 + m_t*(m_t*IT_0002*IT_0006*U_st_00 + mu_h*U_st_10)));
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*conj(V_td)*IT_0155;
    const ccomplex_t IT_0157 = IT_0006*IT_0008;
    const ccomplex_t IT_0158 = m_t*V_td*V_u2*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0159 = IT_0157*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = V_td*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0164 = IT_0008*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0162 + IT_0167;
    const ccomplex_t IT_0169 = IT_0156*IT_0168;
    const ccomplex_t IT_0170 = IT_0154*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = cpow((-2)*s_13 + IT_0001 + IT_0014 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0173 = A_t*IT_0002*U_st_11 + (-1.4142135623731)*e_em
      *IT_0008*(m_W*IT_0002*IT_0004*U_st_01 + (-0.5)*IT_0000*(IT_0001*IT_0003
      *IT_0004*U_st_01 + m_t*(m_t*IT_0002*IT_0006*U_st_01 + mu_h*U_st_11)));
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*conj(V_td)*IT_0173;
    const ccomplex_t IT_0175 = m_t*V_td*V_u2*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0176 = IT_0157*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = (-0.5)*IT_0178;
    const ccomplex_t IT_0180 = V_td*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0181 = IT_0008*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = -IT_0182;
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = IT_0179 + IT_0184;
    const ccomplex_t IT_0186 = IT_0174*IT_0185;
    const ccomplex_t IT_0187 = IT_0172*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0000*IT_0003*IT_0008;
    const ccomplex_t IT_0190 = (-0.5)*IT_0189;
    const ccomplex_t IT_0191 = IT_0126*IT_0190;
    const ccomplex_t IT_0192 = IT_0105*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = m_d*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0000*IT_0003*IT_0008;
    const ccomplex_t IT_0196 = (-0.5)*IT_0195;
    const ccomplex_t IT_0197 = IT_0095*IT_0196;
    const ccomplex_t IT_0198 = IT_0074*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = m_d*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0000*IT_0003*IT_0008;
    const ccomplex_t IT_0202 = (-0.5)*IT_0201;
    const ccomplex_t IT_0203 = IT_0151*IT_0202;
    const ccomplex_t IT_0204 = IT_0130*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = m_d*IT_0205;
    const ccomplex_t IT_0207 = IT_0017 + IT_0025 + IT_0033 + IT_0042 + IT_0072
       + IT_0097 + -IT_0103 + IT_0128 + IT_0153 + -IT_0171 + -IT_0188 + -IT_0194
       + -IT_0200 + -IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0000*IT_0003*IT_0008*conj(U_st_01);
    const ccomplex_t IT_0209 = (-0.5)*IT_0208;
    const ccomplex_t IT_0210 = IT_0174*IT_0209;
    const ccomplex_t IT_0211 = IT_0172*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0000*IT_0003*IT_0008;
    const ccomplex_t IT_0214 = (-0.5)*IT_0213;
    const ccomplex_t IT_0215 = IT_0036*IT_0214;
    const ccomplex_t IT_0216 = IT_0040*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_cd*e_em*IT_0000*IT_0003*IT_0008;
    const ccomplex_t IT_0219 = (-0.5)*IT_0218;
    const ccomplex_t IT_0220 = IT_0010*IT_0219;
    const ccomplex_t IT_0221 = IT_0015*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = IT_0004*IT_0008;
    const ccomplex_t IT_0224 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0225 = IT_0223*IT_0224;
    const ccomplex_t IT_0226 = (-2)*IT_0225;
    const ccomplex_t IT_0227 = IT_0004*IT_0047;
    const ccomplex_t IT_0228 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0229 = IT_0227*IT_0228;
    const ccomplex_t IT_0230 = 1.4142135623731*IT_0229;
    const ccomplex_t IT_0231 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0232 = IT_0223*IT_0231;
    const ccomplex_t IT_0233 = 1.4142135623731*IT_0232;
    const ccomplex_t IT_0234 = IT_0230 + IT_0233;
    const ccomplex_t IT_0235 = IT_0226 + IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = 0.5*IT_0236;
    const ccomplex_t IT_0238 = IT_0190*IT_0237;
    const ccomplex_t IT_0239 = IT_0104*IT_0105*IT_0238;
    const ccomplex_t IT_0240 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0241 = IT_0223*IT_0240;
    const ccomplex_t IT_0242 = (-2)*IT_0241;
    const ccomplex_t IT_0243 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0244 = IT_0227*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0247 = IT_0223*IT_0246;
    const ccomplex_t IT_0248 = 1.4142135623731*IT_0247;
    const ccomplex_t IT_0249 = IT_0245 + IT_0248;
    const ccomplex_t IT_0250 = IT_0242 + IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = 0.5*IT_0251;
    const ccomplex_t IT_0253 = IT_0196*IT_0252;
    const ccomplex_t IT_0254 = IT_0073*IT_0074*IT_0253;
    const ccomplex_t IT_0255 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0256 = IT_0223*IT_0255;
    const ccomplex_t IT_0257 = (-2)*IT_0256;
    const ccomplex_t IT_0258 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0259 = IT_0227*IT_0258;
    const ccomplex_t IT_0260 = 1.4142135623731*IT_0259;
    const ccomplex_t IT_0261 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0262 = IT_0223*IT_0261;
    const ccomplex_t IT_0263 = 1.4142135623731*IT_0262;
    const ccomplex_t IT_0264 = IT_0260 + IT_0263;
    const ccomplex_t IT_0265 = IT_0257 + IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = 0.5*IT_0266;
    const ccomplex_t IT_0268 = IT_0202*IT_0267;
    const ccomplex_t IT_0269 = IT_0129*IT_0130*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0000*IT_0003*IT_0008*conj(U_st_00);
    const ccomplex_t IT_0271 = (-0.5)*IT_0270;
    const ccomplex_t IT_0272 = IT_0156*IT_0271;
    const ccomplex_t IT_0273 = IT_0154*IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*IT_0273;
    const ccomplex_t IT_0275 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0276 = IT_0223*IT_0275;
    const ccomplex_t IT_0277 = 1.4142135623731*IT_0276;
    const ccomplex_t IT_0278 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0279 = IT_0227*IT_0278;
    const ccomplex_t IT_0280 = 1.4142135623731*IT_0279;
    const ccomplex_t IT_0281 = IT_0277 + IT_0280;
    const ccomplex_t IT_0282 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0283 = IT_0223*IT_0282;
    const ccomplex_t IT_0284 = (-2)*IT_0283;
    const ccomplex_t IT_0285 = IT_0281 + IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = 0.5*IT_0286;
    const ccomplex_t IT_0288 = IT_0055*IT_0287;
    const ccomplex_t IT_0289 = IT_0045*IT_0288;
    const ccomplex_t IT_0290 = (0 + _Complex_I*1)*IT_0289;
    const ccomplex_t IT_0291 = m_d*IT_0290;
    const ccomplex_t IT_0292 = IT_0113*IT_0237;
    const ccomplex_t IT_0293 = IT_0105*IT_0292;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*IT_0293;
    const ccomplex_t IT_0295 = m_d*IT_0294;
    const ccomplex_t IT_0296 = IT_0082*IT_0252;
    const ccomplex_t IT_0297 = IT_0074*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = m_d*IT_0298;
    const ccomplex_t IT_0300 = IT_0138*IT_0267;
    const ccomplex_t IT_0301 = IT_0130*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*IT_0301;
    const ccomplex_t IT_0303 = m_d*IT_0302;
    const ccomplex_t IT_0304 = IT_0099*IT_0287;
    const ccomplex_t IT_0305 = IT_0043*IT_0045*IT_0304;
    const ccomplex_t IT_0306 = -IT_0212 + -IT_0217 + -IT_0222 + IT_0239 +
       IT_0254 + IT_0269 + -IT_0274 + -IT_0291 + -IT_0295 + -IT_0299 + -IT_0303 
      + IT_0305;
    const ccomplex_t IT_0307 = m_d*m_C_2;
    const ccomplex_t IT_0308 = -IT_0290 + -IT_0294 + -IT_0298 + -IT_0302;
    const ccomplex_t IT_0309 = s_23*m_C_2;
    const ccomplex_t IT_0310 = -IT_0102 + -IT_0193 + -IT_0199 + -IT_0205;
    const ccomplex_t IT_0311 = m_d*s_12;
    const ccomplex_t IT_0312 = 6*s_13;
    const ccomplex_t IT_0313 = (-6)*IT_0306;
    const ccomplex_t IT_0314 = m_d*m_C_2*IT_0044;
    const ccomplex_t IT_0315 = s_13*IT_0044;
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

