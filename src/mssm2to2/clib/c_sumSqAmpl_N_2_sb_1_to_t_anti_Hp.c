#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_sb_1_to_t_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_sb_1_to_t_anti_Hp(
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
    const creal_t m_sb_1 = param->m_sb_1;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0001 = pow(m_t, 2);
    const ccomplex_t IT_0002 = pow(m_sb_1, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = IT_0004*IT_0006;
    const ccomplex_t IT_0008 = N_d2*e_em*U_Wm2;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = (-2)*IT_0009;
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = N_B2*U_d2*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = N_W2*U_d2*e_em;
    const ccomplex_t IT_0018 = IT_0007*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = IT_0016 + IT_0019;
    const ccomplex_t IT_0021 = IT_0010 + IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = pow(m_W, -1);
    const ccomplex_t IT_0025 = cpow(IT_0004, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0006*IT_0024*IT_0025*U_sb_00;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0023*IT_0027;
    const ccomplex_t IT_0029 = IT_0000*IT_0003*IT_0028;
    const ccomplex_t IT_0030 = pow(m_b, 2);
    const ccomplex_t IT_0031 = pow(m_N_2, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + -IT_0002 + IT_0030 + -IT_0031 
      + -reg_prop, -1);
    const ccomplex_t IT_0033 = cos(beta);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0006*IT_0024*IT_0025*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0033, -1);
    const ccomplex_t IT_0037 = IT_0006*IT_0036;
    const ccomplex_t IT_0038 = m_b*conj(N_d2)*e_em*IT_0024*U_sb_10;
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0042 = IT_0006*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0045 = IT_0012*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0040 + -IT_0043 +
       0.333333333333333*IT_0046);
    const ccomplex_t IT_0048 = 3*IT_0047;
    const ccomplex_t IT_0049 = 0.333333333333333*IT_0048;
    const ccomplex_t IT_0050 = 3*IT_0049;
    const ccomplex_t IT_0051 = (-0.333333333333333)*IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = 3*IT_0052;
    const ccomplex_t IT_0054 = 0.166666666666667*IT_0053;
    const ccomplex_t IT_0055 = IT_0035*IT_0054;
    const ccomplex_t IT_0056 = IT_0032*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = m_N_2*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0060 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0061 = IT_0012*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = m_b*N_d2*e_em*IT_0024*U_sb_00;
    const ccomplex_t IT_0064 = IT_0037*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0062 + 1.5*IT_0065);
    const ccomplex_t IT_0067 = 0.666666666666667*IT_0066;
    const ccomplex_t IT_0068 = 1.5*IT_0067;
    const ccomplex_t IT_0069 = 0.666666666666667*IT_0068;
    const ccomplex_t IT_0070 = 1.5*IT_0069;
    const ccomplex_t IT_0071 = 0.666666666666667*IT_0070;
    const ccomplex_t IT_0072 = 1.5*IT_0071;
    const ccomplex_t IT_0073 = 0.666666666666667*IT_0072;
    const ccomplex_t IT_0074 = 1.5*IT_0073;
    const ccomplex_t IT_0075 = (-0.333333333333333)*IT_0074;
    const ccomplex_t IT_0076 = IT_0035*IT_0075;
    const ccomplex_t IT_0077 = IT_0032*IT_0059*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0079 = N_d2*e_em*U_Wm1;
    const ccomplex_t IT_0080 = IT_0007*IT_0079;
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = N_B2*U_d1*e_em;
    const ccomplex_t IT_0083 = IT_0013*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = N_W2*U_d1*e_em;
    const ccomplex_t IT_0086 = IT_0007*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = IT_0084 + IT_0087;
    const ccomplex_t IT_0089 = IT_0081 + IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = m_b*conj(U_d1)*V_tb*e_em*IT_0024*U_sb_10;
    const ccomplex_t IT_0093 = IT_0037*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0096 = IT_0006*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0094 + (-2)*IT_0096);
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = IT_0091*IT_0098;
    const ccomplex_t IT_0100 = IT_0078*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = m_t*IT_0101;
    const ccomplex_t IT_0103 = m_b*conj(U_d2)*V_tb*e_em*IT_0024*U_sb_10;
    const ccomplex_t IT_0104 = IT_0037*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0107 = IT_0006*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0105 + (-2)*IT_0107);
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = IT_0023*IT_0109;
    const ccomplex_t IT_0111 = IT_0003*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = m_t*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0006*IT_0024*IT_0025*U_sb_00;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = IT_0091*IT_0116;
    const ccomplex_t IT_0118 = IT_0078*IT_0114*IT_0117;
    const ccomplex_t IT_0119 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0120 = IT_0012*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = IT_0006*IT_0025;
    const ccomplex_t IT_0123 = m_t*N_u2*e_em*IT_0024*U_st_00;
    const ccomplex_t IT_0124 = IT_0122*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0121 + (-0.75)*IT_0125);
    const ccomplex_t IT_0127 = 0.666666666666667*IT_0126;
    const ccomplex_t IT_0128 = m_t*U_sb_00*(m_t*IT_0025*IT_0033*conj(U_st_00) 
      + mu_h*conj(U_st_10)) + m_b*(m_b*IT_0004*IT_0036*U_sb_00*conj(U_st_00) +
       U_sb_10*(mu_h*conj(U_st_00) + m_t*IT_0025*IT_0036*conj(U_st_10)));
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0004*U_sb_10
      *conj(U_st_00) + (-1.4142135623731)*e_em*IT_0006*((-0.5)*IT_0024*IT_0128 +
       m_W*IT_0004*IT_0033*U_sb_00*conj(U_st_00)) + A_t*IT_0033*U_sb_00*conj
      (U_st_10));
    const ccomplex_t IT_0130 = IT_0127*IT_0129;
    const ccomplex_t IT_0131 = cpow((-2)*s_13 + IT_0001 + IT_0031 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0132 = IT_0130*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0135 = IT_0012*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = m_t*N_u2*e_em*IT_0024*U_st_01;
    const ccomplex_t IT_0138 = IT_0122*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0136 + (-0.75)*IT_0139);
    const ccomplex_t IT_0141 = 0.666666666666667*IT_0140;
    const ccomplex_t IT_0142 = m_t*U_sb_00*(m_t*IT_0025*IT_0033*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0004*IT_0036*U_sb_00*conj(U_st_01) +
       U_sb_10*(mu_h*conj(U_st_01) + m_t*IT_0025*IT_0036*conj(U_st_11)));
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0004*U_sb_10
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0006*((-0.5)*IT_0024*IT_0142 +
       m_W*IT_0004*IT_0033*U_sb_00*conj(U_st_01)) + A_t*IT_0033*U_sb_00*conj
      (U_st_11));
    const ccomplex_t IT_0144 = IT_0141*IT_0143;
    const ccomplex_t IT_0145 = cpow((-2)*s_13 + IT_0001 + IT_0031 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0146 = IT_0144*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = IT_0029 + IT_0058 + -IT_0077 + -IT_0102 + 
      -IT_0113 + IT_0118 + IT_0133 + IT_0147;
    const ccomplex_t IT_0149 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0150 = IT_0012*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0153 = IT_0006*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = m_t*conj(N_u2)*e_em*IT_0024*U_st_11;
    const ccomplex_t IT_0156 = IT_0122*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0151 + 3*IT_0154 + 3
      *IT_0157);
    const ccomplex_t IT_0159 = 0.166666666666667*IT_0158;
    const ccomplex_t IT_0160 = IT_0143*IT_0159;
    const ccomplex_t IT_0161 = IT_0145*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = IT_0006*IT_0033;
    const ccomplex_t IT_0164 = conj(N_u2)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0165 = IT_0163*IT_0164;
    const ccomplex_t IT_0166 = 2*IT_0165;
    const ccomplex_t IT_0167 = IT_0012*IT_0033;
    const ccomplex_t IT_0168 = conj(N_B2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0169 = IT_0167*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = conj(N_W2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0172 = IT_0163*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = IT_0170 + IT_0173;
    const ccomplex_t IT_0175 = IT_0166 + IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = 0.5*IT_0176;
    const ccomplex_t IT_0178 = IT_0098*IT_0177;
    const ccomplex_t IT_0179 = IT_0078*IT_0114*IT_0178;
    const ccomplex_t IT_0180 = conj(N_u2)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0181 = IT_0163*IT_0180;
    const ccomplex_t IT_0182 = 2*IT_0181;
    const ccomplex_t IT_0183 = conj(N_B2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0184 = IT_0167*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = conj(N_W2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0187 = IT_0163*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = IT_0185 + IT_0188;
    const ccomplex_t IT_0190 = IT_0182 + IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = IT_0109*IT_0192;
    const ccomplex_t IT_0194 = IT_0000*IT_0003*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0004*IT_0006*IT_0024*IT_0036;
    const ccomplex_t IT_0196 = 0.5*IT_0195;
    const ccomplex_t IT_0197 = IT_0054*IT_0196;
    const ccomplex_t IT_0198 = IT_0032*IT_0059*IT_0197;
    const ccomplex_t IT_0199 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0200 = IT_0012*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0203 = IT_0006*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = m_t*conj(N_u2)*e_em*IT_0024*U_st_10;
    const ccomplex_t IT_0206 = IT_0122*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0201 + 3*IT_0204 + 3
      *IT_0207);
    const ccomplex_t IT_0209 = 0.166666666666667*IT_0208;
    const ccomplex_t IT_0210 = IT_0129*IT_0209;
    const ccomplex_t IT_0211 = IT_0131*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = IT_0116*IT_0177;
    const ccomplex_t IT_0214 = IT_0078*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = m_t*IT_0215;
    const ccomplex_t IT_0217 = IT_0027*IT_0192;
    const ccomplex_t IT_0218 = IT_0003*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = m_t*IT_0219;
    const ccomplex_t IT_0221 = IT_0075*IT_0196;
    const ccomplex_t IT_0222 = IT_0032*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = m_N_2*IT_0223;
    const ccomplex_t IT_0225 = -IT_0162 + IT_0179 + IT_0194 + IT_0198 + 
      -IT_0212 + -IT_0216 + -IT_0220 + -IT_0224;
    const ccomplex_t IT_0226 = 6*s_13;
    const ccomplex_t IT_0227 = m_t*m_N_2;
    const ccomplex_t IT_0228 = 6*IT_0148;
    const ccomplex_t IT_0229 = IT_0101 + IT_0112 + -IT_0223;
    const ccomplex_t IT_0230 = m_t*s_12;
    const ccomplex_t IT_0231 = IT_0057 + IT_0215 + IT_0219;
    const ccomplex_t IT_0232 = s_23*m_N_2;
    const ccomplex_t IT_0233 = 6*conj(IT_0148);
    const ccomplex_t IT_0234 = 6*conj(IT_0225);
    const ccomplex_t IT_0235 = m_t*m_N_2*IT_0002;
    const ccomplex_t IT_0236 = s_13*IT_0002;
    const ccomplex_t IT_0237 = (-6)*IT_0236;
    const ccomplex_t IT_0238 = s_12*s_23;
    const ccomplex_t IT_0239 = 12*IT_0238;
    const ccomplex_t IT_0240 = IT_0237 + IT_0239;
    const ccomplex_t IT_0241 = 6*IT_0225;
    const ccomplex_t IT_0242 = 0.166666666666667*IT_0230;
    const ccomplex_t IT_0243 = 0.166666666666667*IT_0232;
    const ccomplex_t IT_0244 = (IT_0148*conj(IT_0148) + IT_0225*conj(IT_0225))
      *IT_0226 + conj(IT_0225)*IT_0227*IT_0228 + (IT_0225*IT_0227 + IT_0229
      *IT_0230 + IT_0231*IT_0232)*IT_0233 + (IT_0230*IT_0231 + IT_0229*IT_0232)
      *IT_0234 + 6*conj(IT_0231)*(IT_0229*IT_0235 + 0.166666666666667*IT_0231
      *IT_0240 + IT_0241*IT_0242 + IT_0228*IT_0243) + 6*conj(IT_0229)*(IT_0231
      *IT_0235 + 0.166666666666667*IT_0229*IT_0240 + IT_0228*IT_0242 + IT_0241
      *IT_0243);
    return create_ccomplex_return(IT_0244);
}

