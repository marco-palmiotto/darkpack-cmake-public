#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_sb_2_to_t_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_sb_2_to_t_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
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
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = m_t*N_u2*e_em*IT_0010*U_st_00;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0004 + (-0.75)*IT_0013);
    const ccomplex_t IT_0015 = 0.666666666666667*IT_0014;
    const ccomplex_t IT_0016 = cos(beta);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = m_t*U_sb_01*(m_t*IT_0006*IT_0016*conj(U_st_00) 
      + mu_h*conj(U_st_10)) + m_b*(m_b*IT_0005*IT_0017*U_sb_01*conj(U_st_00) +
       U_sb_11*(mu_h*conj(U_st_00) + m_t*IT_0006*IT_0017*conj(U_st_10)));
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*U_sb_11
      *conj(U_st_00) + (-1.4142135623731)*e_em*IT_0008*((-0.5)*IT_0010*IT_0018 +
       m_W*IT_0005*IT_0016*U_sb_01*conj(U_st_00)) + A_t*IT_0016*U_sb_01*conj
      (U_st_10));
    const ccomplex_t IT_0020 = IT_0015*IT_0019;
    const ccomplex_t IT_0021 = pow(m_t, 2);
    const ccomplex_t IT_0022 = pow(m_N_2, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0021 + IT_0022 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = m_t*N_u2*e_em*IT_0010*U_st_01;
    const ccomplex_t IT_0030 = IT_0009*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + (-0.75)*IT_0031);
    const ccomplex_t IT_0033 = 0.666666666666667*IT_0032;
    const ccomplex_t IT_0034 = m_t*U_sb_01*(m_t*IT_0006*IT_0016*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0005*IT_0017*U_sb_01*conj(U_st_01) +
       U_sb_11*(mu_h*conj(U_st_01) + m_t*IT_0006*IT_0017*conj(U_st_11)));
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*U_sb_11
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0008*((-0.5)*IT_0010*IT_0034 +
       m_W*IT_0005*IT_0016*U_sb_01*conj(U_st_01)) + A_t*IT_0016*U_sb_01*conj
      (U_st_11));
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0021 + IT_0022 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = pow(m_sb_2, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0021 + IT_0040 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0005*IT_0008;
    const ccomplex_t IT_0043 = N_d2*e_em*U_Wm1;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (-2)*IT_0044;
    const ccomplex_t IT_0046 = IT_0001*IT_0005;
    const ccomplex_t IT_0047 = N_B2*U_d1*e_em;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = N_W2*U_d1*e_em;
    const ccomplex_t IT_0051 = IT_0042*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = IT_0049 + IT_0052;
    const ccomplex_t IT_0054 = IT_0045 + IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = IT_0008*IT_0017;
    const ccomplex_t IT_0058 = m_b*conj(U_d1)*V_tb*e_em*IT_0010*U_sb_11;
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0062 = IT_0008*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0060 + (-2)*IT_0062);
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = IT_0056*IT_0064;
    const ccomplex_t IT_0066 = IT_0041*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = m_t*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0021 + IT_0040 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0070 = N_d2*e_em*U_Wm2;
    const ccomplex_t IT_0071 = IT_0042*IT_0070;
    const ccomplex_t IT_0072 = (-2)*IT_0071;
    const ccomplex_t IT_0073 = N_B2*U_d2*e_em;
    const ccomplex_t IT_0074 = IT_0046*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = N_W2*U_d2*e_em;
    const ccomplex_t IT_0077 = IT_0042*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = IT_0075 + IT_0078;
    const ccomplex_t IT_0080 = IT_0072 + IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = m_b*conj(U_d2)*V_tb*e_em*IT_0010*U_sb_11;
    const ccomplex_t IT_0084 = IT_0057*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0087 = IT_0008*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0085 + (-2)*IT_0087);
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0082*IT_0089;
    const ccomplex_t IT_0091 = IT_0069*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = m_t*IT_0092;
    const ccomplex_t IT_0094 = pow(m_b, 2);
    const ccomplex_t IT_0095 = cpow((-2)*s_12 + -IT_0022 + -IT_0040 + IT_0094 
      + -reg_prop, -1);
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0006*IT_0008*IT_0010*IT_0016;
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = m_b*conj(N_d2)*e_em*IT_0010*U_sb_11;
    const ccomplex_t IT_0099 = IT_0057*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0102 = IT_0008*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0105 = IT_0001*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0100 + -IT_0103 +
       0.333333333333333*IT_0106);
    const ccomplex_t IT_0108 = 3*IT_0107;
    const ccomplex_t IT_0109 = 0.166666666666667*IT_0108;
    const ccomplex_t IT_0110 = IT_0097*IT_0109;
    const ccomplex_t IT_0111 = IT_0095*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = m_N_2*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0006*IT_0008*IT_0010*U_sb_01;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = IT_0056*IT_0116;
    const ccomplex_t IT_0118 = IT_0041*IT_0114*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0006*IT_0008*IT_0010*U_sb_01;
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = IT_0082*IT_0121;
    const ccomplex_t IT_0123 = IT_0069*IT_0119*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0125 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0126 = IT_0001*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = m_b*N_d2*e_em*IT_0010*U_sb_01;
    const ccomplex_t IT_0129 = IT_0057*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + 1.5*IT_0130);
    const ccomplex_t IT_0132 = (-0.333333333333333)*IT_0131;
    const ccomplex_t IT_0133 = IT_0097*IT_0132;
    const ccomplex_t IT_0134 = IT_0095*IT_0124*IT_0133;
    const ccomplex_t IT_0135 = IT_0025 + IT_0039 + -IT_0068 + -IT_0093 +
       IT_0113 + IT_0118 + IT_0123 + -IT_0134;
    const ccomplex_t IT_0136 = IT_0008*IT_0016;
    const ccomplex_t IT_0137 = conj(N_u2)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0138 = IT_0136*IT_0137;
    const ccomplex_t IT_0139 = 2*IT_0138;
    const ccomplex_t IT_0140 = IT_0001*IT_0016;
    const ccomplex_t IT_0141 = conj(N_B2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0142 = IT_0140*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = conj(N_W2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0145 = IT_0136*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = IT_0143 + IT_0146;
    const ccomplex_t IT_0148 = IT_0139 + IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = IT_0064*IT_0150;
    const ccomplex_t IT_0152 = IT_0041*IT_0114*IT_0151;
    const ccomplex_t IT_0153 = conj(N_u2)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0154 = IT_0136*IT_0153;
    const ccomplex_t IT_0155 = 2*IT_0154;
    const ccomplex_t IT_0156 = conj(N_B2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0157 = IT_0140*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = conj(N_W2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0160 = IT_0136*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = IT_0158 + IT_0161;
    const ccomplex_t IT_0163 = IT_0155 + IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = 0.5*IT_0164;
    const ccomplex_t IT_0166 = IT_0089*IT_0165;
    const ccomplex_t IT_0167 = IT_0069*IT_0119*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0169 = IT_0001*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0172 = IT_0008*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = m_t*conj(N_u2)*e_em*IT_0010*U_st_10;
    const ccomplex_t IT_0175 = IT_0009*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*(IT_0170 + 3*IT_0173 + 3
      *IT_0176);
    const ccomplex_t IT_0178 = 0.166666666666667*IT_0177;
    const ccomplex_t IT_0179 = IT_0019*IT_0178;
    const ccomplex_t IT_0180 = IT_0023*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0183 = IT_0001*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0186 = IT_0008*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = m_t*conj(N_u2)*e_em*IT_0010*U_st_11;
    const ccomplex_t IT_0189 = IT_0009*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*(IT_0184 + 3*IT_0187 + 3
      *IT_0190);
    const ccomplex_t IT_0192 = 0.166666666666667*IT_0191;
    const ccomplex_t IT_0193 = IT_0035*IT_0192;
    const ccomplex_t IT_0194 = IT_0037*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0005*IT_0008*IT_0010*IT_0017;
    const ccomplex_t IT_0197 = 0.5*IT_0196;
    const ccomplex_t IT_0198 = IT_0109*IT_0197;
    const ccomplex_t IT_0199 = IT_0095*IT_0124*IT_0198;
    const ccomplex_t IT_0200 = IT_0121*IT_0165;
    const ccomplex_t IT_0201 = IT_0069*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = m_t*IT_0202;
    const ccomplex_t IT_0204 = IT_0116*IT_0150;
    const ccomplex_t IT_0205 = IT_0041*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = m_t*IT_0206;
    const ccomplex_t IT_0208 = IT_0132*IT_0197;
    const ccomplex_t IT_0209 = IT_0095*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = m_N_2*IT_0210;
    const ccomplex_t IT_0212 = IT_0152 + IT_0167 + -IT_0181 + -IT_0195 +
       IT_0199 + -IT_0203 + -IT_0207 + -IT_0211;
    const ccomplex_t IT_0213 = 6*s_13;
    const ccomplex_t IT_0214 = m_t*m_N_2;
    const ccomplex_t IT_0215 = 6*IT_0135;
    const ccomplex_t IT_0216 = IT_0067 + IT_0092 + -IT_0210;
    const ccomplex_t IT_0217 = m_t*s_12;
    const ccomplex_t IT_0218 = IT_0112 + IT_0202 + IT_0206;
    const ccomplex_t IT_0219 = s_23*m_N_2;
    const ccomplex_t IT_0220 = 6*conj(IT_0135);
    const ccomplex_t IT_0221 = 6*conj(IT_0212);
    const ccomplex_t IT_0222 = m_t*m_N_2*IT_0040;
    const ccomplex_t IT_0223 = s_13*IT_0040;
    const ccomplex_t IT_0224 = (-6)*IT_0223;
    const ccomplex_t IT_0225 = s_12*s_23;
    const ccomplex_t IT_0226 = 12*IT_0225;
    const ccomplex_t IT_0227 = IT_0224 + IT_0226;
    const ccomplex_t IT_0228 = 6*IT_0212;
    const ccomplex_t IT_0229 = 0.166666666666667*IT_0217;
    const ccomplex_t IT_0230 = 0.166666666666667*IT_0219;
    const ccomplex_t IT_0231 = (IT_0135*conj(IT_0135) + IT_0212*conj(IT_0212))
      *IT_0213 + conj(IT_0212)*IT_0214*IT_0215 + (IT_0212*IT_0214 + IT_0216
      *IT_0217 + IT_0218*IT_0219)*IT_0220 + (IT_0217*IT_0218 + IT_0216*IT_0219)
      *IT_0221 + 6*conj(IT_0218)*(IT_0216*IT_0222 + 0.166666666666667*IT_0218
      *IT_0227 + IT_0228*IT_0229 + IT_0215*IT_0230) + 6*conj(IT_0216)*(IT_0218
      *IT_0222 + 0.166666666666667*IT_0216*IT_0227 + IT_0215*IT_0229 + IT_0228
      *IT_0230);
    return create_ccomplex_return(IT_0231);
}

