#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_sb_2_to_t_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_sb_2_to_t_anti_Hp(
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
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
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
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_N_1, 2);
    const ccomplex_t IT_0002 = pow(m_sb_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0005*IT_0006*IT_0008*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = cos(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0006, -1);
    const ccomplex_t IT_0019 = IT_0010*IT_0018;
    const ccomplex_t IT_0020 = m_b*N_d1*e_em*IT_0005*U_sb_01;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0017 + 1.5*IT_0022);
    const ccomplex_t IT_0024 = (-0.333333333333333)*IT_0023;
    const ccomplex_t IT_0025 = IT_0012*IT_0024;
    const ccomplex_t IT_0026 = IT_0003*IT_0004*IT_0025;
    const ccomplex_t IT_0027 = pow(m_t, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0002 + IT_0027 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0007*IT_0010;
    const ccomplex_t IT_0030 = N_d1*e_em*U_Wm2;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (-2)*IT_0031;
    const ccomplex_t IT_0033 = IT_0007*IT_0014;
    const ccomplex_t IT_0034 = N_B1*U_d2*e_em;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_W1*U_d2*e_em;
    const ccomplex_t IT_0038 = IT_0029*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = IT_0036 + IT_0039;
    const ccomplex_t IT_0041 = IT_0032 + IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = m_b*conj(U_d2)*V_tb*e_em*IT_0005*U_sb_11;
    const ccomplex_t IT_0045 = IT_0019*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0048 = IT_0010*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0046 + (-2)*IT_0048);
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0043*IT_0050;
    const ccomplex_t IT_0052 = IT_0028*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = m_t*IT_0053;
    const ccomplex_t IT_0055 = m_b*conj(N_d1)*e_em*IT_0005*U_sb_11;
    const ccomplex_t IT_0056 = IT_0019*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0059 = IT_0010*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0062 = IT_0014*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0057 + -IT_0060 +
       0.333333333333333*IT_0063);
    const ccomplex_t IT_0065 = 3*IT_0064;
    const ccomplex_t IT_0066 = 0.166666666666667*IT_0065;
    const ccomplex_t IT_0067 = IT_0012*IT_0066;
    const ccomplex_t IT_0068 = IT_0003*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = m_N_1*IT_0069;
    const ccomplex_t IT_0071 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0072 = IT_0014*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = IT_0008*IT_0010;
    const ccomplex_t IT_0075 = m_t*N_u1*e_em*IT_0005*U_st_00;
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0073 + (-0.75)*IT_0077);
    const ccomplex_t IT_0079 = 0.666666666666667*IT_0078;
    const ccomplex_t IT_0080 = m_t*U_sb_01*(m_t*IT_0006*IT_0008*conj(U_st_00) 
      + mu_h*conj(U_st_10)) + m_b*(m_b*IT_0007*IT_0018*U_sb_01*conj(U_st_00) +
       U_sb_11*(mu_h*conj(U_st_00) + m_t*IT_0008*IT_0018*conj(U_st_10)));
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0007*U_sb_11
      *conj(U_st_00) + (-1.4142135623731)*e_em*IT_0010*((-0.5)*IT_0005*IT_0080 +
       m_W*IT_0006*IT_0007*U_sb_01*conj(U_st_00)) + A_t*IT_0006*U_sb_01*conj
      (U_st_10));
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_13 + IT_0001 + IT_0027 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0087 = IT_0014*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = m_t*N_u1*e_em*IT_0005*U_st_01;
    const ccomplex_t IT_0090 = IT_0074*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + (-0.75)*IT_0091);
    const ccomplex_t IT_0093 = 0.666666666666667*IT_0092;
    const ccomplex_t IT_0094 = m_t*U_sb_01*(m_t*IT_0006*IT_0008*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0007*IT_0018*U_sb_01*conj(U_st_01) +
       U_sb_11*(mu_h*conj(U_st_01) + m_t*IT_0008*IT_0018*conj(U_st_11)));
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0007*U_sb_11
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0010*((-0.5)*IT_0005*IT_0094 +
       m_W*IT_0006*IT_0007*U_sb_01*conj(U_st_01)) + A_t*IT_0006*U_sb_01*conj
      (U_st_11));
    const ccomplex_t IT_0096 = IT_0093*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_13 + IT_0001 + IT_0027 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = cpow((-2)*s_23 + IT_0002 + IT_0027 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0101 = N_d1*e_em*U_Wm1;
    const ccomplex_t IT_0102 = IT_0029*IT_0101;
    const ccomplex_t IT_0103 = (-2)*IT_0102;
    const ccomplex_t IT_0104 = N_B1*U_d1*e_em;
    const ccomplex_t IT_0105 = IT_0033*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = N_W1*U_d1*e_em;
    const ccomplex_t IT_0108 = IT_0029*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = IT_0106 + IT_0109;
    const ccomplex_t IT_0111 = IT_0103 + IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = m_b*conj(U_d1)*V_tb*e_em*IT_0005*U_sb_11;
    const ccomplex_t IT_0115 = IT_0019*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0118 = IT_0010*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0116 + (-2)*IT_0118);
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = IT_0113*IT_0120;
    const ccomplex_t IT_0122 = IT_0100*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = m_t*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0005*IT_0008*IT_0010*U_sb_01;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = IT_0113*IT_0127;
    const ccomplex_t IT_0129 = IT_0100*IT_0125*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0005*IT_0008*IT_0010*U_sb_01;
    const ccomplex_t IT_0132 = 0.5*IT_0131;
    const ccomplex_t IT_0133 = IT_0043*IT_0132;
    const ccomplex_t IT_0134 = IT_0028*IT_0130*IT_0133;
    const ccomplex_t IT_0135 = -IT_0026 + -IT_0054 + IT_0070 + IT_0085 +
       IT_0099 + -IT_0124 + IT_0129 + IT_0134;
    const ccomplex_t IT_0136 = IT_0006*IT_0010;
    const ccomplex_t IT_0137 = conj(N_u1)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0138 = IT_0136*IT_0137;
    const ccomplex_t IT_0139 = 2*IT_0138;
    const ccomplex_t IT_0140 = IT_0006*IT_0014;
    const ccomplex_t IT_0141 = conj(N_B1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0142 = IT_0140*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = conj(N_W1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0145 = IT_0136*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = IT_0143 + IT_0146;
    const ccomplex_t IT_0148 = IT_0139 + IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = IT_0120*IT_0150;
    const ccomplex_t IT_0152 = IT_0100*IT_0125*IT_0151;
    const ccomplex_t IT_0153 = conj(N_u1)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0154 = IT_0136*IT_0153;
    const ccomplex_t IT_0155 = 2*IT_0154;
    const ccomplex_t IT_0156 = conj(N_B1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0157 = IT_0140*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = conj(N_W1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0160 = IT_0136*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = IT_0158 + IT_0161;
    const ccomplex_t IT_0163 = IT_0155 + IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = 0.5*IT_0164;
    const ccomplex_t IT_0166 = IT_0050*IT_0165;
    const ccomplex_t IT_0167 = IT_0028*IT_0130*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0005*IT_0007*IT_0010*IT_0018;
    const ccomplex_t IT_0169 = 0.5*IT_0168;
    const ccomplex_t IT_0170 = IT_0066*IT_0169;
    const ccomplex_t IT_0171 = IT_0003*IT_0004*IT_0170;
    const ccomplex_t IT_0172 = IT_0024*IT_0169;
    const ccomplex_t IT_0173 = IT_0003*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = m_N_1*IT_0174;
    const ccomplex_t IT_0176 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0177 = IT_0014*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0180 = IT_0010*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = m_t*conj(N_u1)*e_em*IT_0005*U_st_10;
    const ccomplex_t IT_0183 = IT_0074*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0178 + 3*IT_0181 + 3
      *IT_0184);
    const ccomplex_t IT_0186 = 0.166666666666667*IT_0185;
    const ccomplex_t IT_0187 = IT_0081*IT_0186;
    const ccomplex_t IT_0188 = IT_0083*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0191 = IT_0014*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0194 = IT_0010*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = m_t*conj(N_u1)*e_em*IT_0005*U_st_11;
    const ccomplex_t IT_0197 = IT_0074*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0192 + 3*IT_0195 + 3
      *IT_0198);
    const ccomplex_t IT_0200 = 0.166666666666667*IT_0199;
    const ccomplex_t IT_0201 = IT_0095*IT_0200;
    const ccomplex_t IT_0202 = IT_0097*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = IT_0127*IT_0150;
    const ccomplex_t IT_0205 = IT_0100*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = m_t*IT_0206;
    const ccomplex_t IT_0208 = IT_0132*IT_0165;
    const ccomplex_t IT_0209 = IT_0028*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = m_t*IT_0210;
    const ccomplex_t IT_0212 = IT_0152 + IT_0167 + IT_0171 + -IT_0175 + 
      -IT_0189 + -IT_0203 + -IT_0207 + -IT_0211;
    const ccomplex_t IT_0213 = 6*s_13;
    const ccomplex_t IT_0214 = m_t*m_N_1;
    const ccomplex_t IT_0215 = 6*IT_0135;
    const ccomplex_t IT_0216 = IT_0053 + IT_0123 + -IT_0174;
    const ccomplex_t IT_0217 = m_t*s_12;
    const ccomplex_t IT_0218 = IT_0069 + IT_0206 + IT_0210;
    const ccomplex_t IT_0219 = s_23*m_N_1;
    const ccomplex_t IT_0220 = 6*conj(IT_0135);
    const ccomplex_t IT_0221 = 6*conj(IT_0212);
    const ccomplex_t IT_0222 = m_t*m_N_1*IT_0002;
    const ccomplex_t IT_0223 = s_13*IT_0002;
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

