#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_sb_2_to_t_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_sb_2_to_t_anti_Hp(
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
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
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
    const ccomplex_t IT_0000 = pow(m_t, 2);
    const ccomplex_t IT_0001 = pow(m_sb_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = N_d3*e_em*U_Wm1;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = (-2)*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0003*IT_0011;
    const ccomplex_t IT_0013 = N_B3*U_d1*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = N_W3*U_d1*e_em;
    const ccomplex_t IT_0017 = IT_0006*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = IT_0015 + IT_0018;
    const ccomplex_t IT_0020 = IT_0009 + IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = cos(beta);
    const ccomplex_t IT_0024 = cpow(IT_0023, -1);
    const ccomplex_t IT_0025 = IT_0005*IT_0024;
    const ccomplex_t IT_0026 = pow(m_W, -1);
    const ccomplex_t IT_0027 = m_b*conj(U_d1)*V_tb*e_em*IT_0026*U_sb_11;
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0031 = IT_0005*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0029 + (-2)*IT_0031);
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = IT_0022*IT_0033;
    const ccomplex_t IT_0035 = IT_0002*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = m_t*IT_0036;
    const ccomplex_t IT_0038 = pow(m_b, 2);
    const ccomplex_t IT_0039 = pow(m_N_3, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_12 + -IT_0001 + IT_0038 + -IT_0039 
      + -reg_prop, -1);
    const ccomplex_t IT_0041 = cpow(IT_0003, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0005*IT_0023*IT_0026*IT_0041;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0045 = IT_0011*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0048 = IT_0005*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = m_b*conj(N_d3)*e_em*IT_0026*U_sb_11;
    const ccomplex_t IT_0051 = IT_0025*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0046 + (-3)*IT_0049 + 3
      *IT_0052);
    const ccomplex_t IT_0054 = 0.166666666666667*IT_0053;
    const ccomplex_t IT_0055 = IT_0043*IT_0054;
    const ccomplex_t IT_0056 = IT_0040*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = m_N_3*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0060 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0061 = N_W3*U_d2*e_em;
    const ccomplex_t IT_0062 = IT_0006*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = N_B3*U_d2*e_em;
    const ccomplex_t IT_0065 = IT_0012*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = IT_0063 + IT_0066;
    const ccomplex_t IT_0068 = N_d3*e_em*U_Wm2;
    const ccomplex_t IT_0069 = IT_0006*IT_0068;
    const ccomplex_t IT_0070 = (-2)*IT_0069;
    const ccomplex_t IT_0071 = IT_0067 + IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0005*IT_0026*IT_0041*U_sb_01;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = IT_0073*IT_0075;
    const ccomplex_t IT_0077 = IT_0059*IT_0060*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_13 + IT_0000 + IT_0039 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0079 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0080 = IT_0011*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = IT_0005*IT_0041;
    const ccomplex_t IT_0083 = m_t*N_u3*e_em*IT_0026*U_st_00;
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0081 + (-0.75)*IT_0085);
    const ccomplex_t IT_0087 = 0.666666666666667*IT_0086;
    const ccomplex_t IT_0088 = m_t*U_sb_01*(m_t*IT_0023*IT_0041*conj(U_st_00) 
      + mu_h*conj(U_st_10)) + m_b*(m_b*IT_0003*IT_0024*U_sb_01*conj(U_st_00) +
       U_sb_11*(mu_h*conj(U_st_00) + m_t*IT_0024*IT_0041*conj(U_st_10)));
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0003*U_sb_11
      *conj(U_st_00) + (-1.4142135623731)*e_em*IT_0005*((-0.5)*IT_0026*IT_0088 +
       m_W*IT_0003*IT_0023*U_sb_01*conj(U_st_00)) + A_t*IT_0023*U_sb_01*conj
      (U_st_10));
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = IT_0078*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = cpow((-2)*s_13 + IT_0000 + IT_0039 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0094 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0095 = IT_0011*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = m_t*N_u3*e_em*IT_0026*U_st_01;
    const ccomplex_t IT_0098 = IT_0082*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + (-0.75)*IT_0099);
    const ccomplex_t IT_0101 = 0.666666666666667*IT_0100;
    const ccomplex_t IT_0102 = m_t*U_sb_01*(m_t*IT_0023*IT_0041*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0003*IT_0024*U_sb_01*conj(U_st_01) +
       U_sb_11*(mu_h*conj(U_st_01) + m_t*IT_0024*IT_0041*conj(U_st_11)));
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0003*U_sb_11
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0005*((-0.5)*IT_0026*IT_0102 +
       m_W*IT_0003*IT_0023*U_sb_01*conj(U_st_01)) + A_t*IT_0023*U_sb_01*conj
      (U_st_11));
    const ccomplex_t IT_0104 = IT_0101*IT_0103;
    const ccomplex_t IT_0105 = IT_0093*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0108 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0109 = IT_0011*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = m_b*N_d3*e_em*IT_0026*U_sb_01;
    const ccomplex_t IT_0112 = IT_0025*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + 1.5*IT_0113);
    const ccomplex_t IT_0115 = (-0.333333333333333)*IT_0114;
    const ccomplex_t IT_0116 = IT_0043*IT_0115;
    const ccomplex_t IT_0117 = IT_0040*IT_0107*IT_0116;
    const ccomplex_t IT_0118 = m_b*conj(U_d2)*V_tb*e_em*IT_0026*U_sb_11;
    const ccomplex_t IT_0119 = IT_0025*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0122 = IT_0005*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0120 + (-2)*IT_0122);
    const ccomplex_t IT_0124 = (-0.5)*IT_0123;
    const ccomplex_t IT_0125 = IT_0073*IT_0124;
    const ccomplex_t IT_0126 = IT_0060*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = m_t*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0005*IT_0026*IT_0041*U_sb_01;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = IT_0022*IT_0131;
    const ccomplex_t IT_0133 = IT_0002*IT_0129*IT_0132;
    const ccomplex_t IT_0134 = -IT_0037 + IT_0058 + IT_0077 + IT_0092 +
       IT_0106 + -IT_0117 + -IT_0128 + IT_0133;
    const ccomplex_t IT_0135 = IT_0005*IT_0023;
    const ccomplex_t IT_0136 = conj(N_u3)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0137 = IT_0135*IT_0136;
    const ccomplex_t IT_0138 = 2*IT_0137;
    const ccomplex_t IT_0139 = IT_0011*IT_0023;
    const ccomplex_t IT_0140 = conj(N_B3)*conj(V_u1)*e_em;
    const ccomplex_t IT_0141 = IT_0139*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = conj(N_W3)*conj(V_u1)*e_em;
    const ccomplex_t IT_0144 = IT_0135*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = IT_0142 + IT_0145;
    const ccomplex_t IT_0147 = IT_0138 + IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = IT_0033*IT_0149;
    const ccomplex_t IT_0151 = IT_0002*IT_0129*IT_0150;
    const ccomplex_t IT_0152 = conj(N_u3)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0153 = IT_0135*IT_0152;
    const ccomplex_t IT_0154 = 2*IT_0153;
    const ccomplex_t IT_0155 = conj(N_B3)*conj(V_u2)*e_em;
    const ccomplex_t IT_0156 = IT_0139*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = conj(N_W3)*conj(V_u2)*e_em;
    const ccomplex_t IT_0159 = IT_0135*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = IT_0157 + IT_0160;
    const ccomplex_t IT_0162 = IT_0154 + IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = IT_0124*IT_0164;
    const ccomplex_t IT_0166 = IT_0059*IT_0060*IT_0165;
    const ccomplex_t IT_0167 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0168 = IT_0011*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0171 = IT_0005*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = m_t*conj(N_u3)*e_em*IT_0026*U_st_11;
    const ccomplex_t IT_0174 = IT_0082*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0169 + 3*IT_0172 + 3
      *IT_0175);
    const ccomplex_t IT_0177 = 0.166666666666667*IT_0176;
    const ccomplex_t IT_0178 = IT_0103*IT_0177;
    const ccomplex_t IT_0179 = IT_0093*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0003*IT_0005*IT_0024*IT_0026;
    const ccomplex_t IT_0182 = 0.5*IT_0181;
    const ccomplex_t IT_0183 = IT_0054*IT_0182;
    const ccomplex_t IT_0184 = IT_0040*IT_0107*IT_0183;
    const ccomplex_t IT_0185 = IT_0115*IT_0182;
    const ccomplex_t IT_0186 = IT_0040*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = m_N_3*IT_0187;
    const ccomplex_t IT_0189 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0190 = IT_0011*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0193 = IT_0005*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = m_t*conj(N_u3)*e_em*IT_0026*U_st_10;
    const ccomplex_t IT_0196 = IT_0082*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*(IT_0191 + 3*IT_0194 + 3
      *IT_0197);
    const ccomplex_t IT_0199 = 0.166666666666667*IT_0198;
    const ccomplex_t IT_0200 = IT_0089*IT_0199;
    const ccomplex_t IT_0201 = IT_0078*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = IT_0131*IT_0149;
    const ccomplex_t IT_0204 = IT_0002*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = m_t*IT_0205;
    const ccomplex_t IT_0207 = IT_0075*IT_0164;
    const ccomplex_t IT_0208 = IT_0060*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = m_t*IT_0209;
    const ccomplex_t IT_0211 = IT_0151 + IT_0166 + -IT_0180 + IT_0184 + 
      -IT_0188 + -IT_0202 + -IT_0206 + -IT_0210;
    const ccomplex_t IT_0212 = 6*s_13;
    const ccomplex_t IT_0213 = m_t*m_N_3;
    const ccomplex_t IT_0214 = 6*IT_0134;
    const ccomplex_t IT_0215 = IT_0036 + IT_0127 + -IT_0187;
    const ccomplex_t IT_0216 = m_t*s_12;
    const ccomplex_t IT_0217 = IT_0057 + IT_0205 + IT_0209;
    const ccomplex_t IT_0218 = s_23*m_N_3;
    const ccomplex_t IT_0219 = 6*conj(IT_0134);
    const ccomplex_t IT_0220 = 6*conj(IT_0211);
    const ccomplex_t IT_0221 = m_t*m_N_3*IT_0001;
    const ccomplex_t IT_0222 = s_13*IT_0001;
    const ccomplex_t IT_0223 = (-6)*IT_0222;
    const ccomplex_t IT_0224 = s_12*s_23;
    const ccomplex_t IT_0225 = 12*IT_0224;
    const ccomplex_t IT_0226 = IT_0223 + IT_0225;
    const ccomplex_t IT_0227 = 6*IT_0211;
    const ccomplex_t IT_0228 = 0.166666666666667*IT_0216;
    const ccomplex_t IT_0229 = 0.166666666666667*IT_0218;
    const ccomplex_t IT_0230 = (IT_0134*conj(IT_0134) + IT_0211*conj(IT_0211))
      *IT_0212 + conj(IT_0211)*IT_0213*IT_0214 + (IT_0211*IT_0213 + IT_0215
      *IT_0216 + IT_0217*IT_0218)*IT_0219 + (IT_0216*IT_0217 + IT_0215*IT_0218)
      *IT_0220 + 6*conj(IT_0217)*(IT_0215*IT_0221 + 0.166666666666667*IT_0217
      *IT_0226 + IT_0227*IT_0228 + IT_0214*IT_0229) + 6*conj(IT_0215)*(IT_0217
      *IT_0221 + 0.166666666666667*IT_0215*IT_0226 + IT_0214*IT_0228 + IT_0227
      *IT_0229);
    return create_ccomplex_return(IT_0230);
}

