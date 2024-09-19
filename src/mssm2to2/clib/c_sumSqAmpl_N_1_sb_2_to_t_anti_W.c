#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_sb_2_to_t_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_sb_2_to_t_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
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
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*U_sb_01*conj(U_st_00);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0001*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = m_t*N_u1*e_em*IT_0012*U_st_00;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0008 + (-0.75)*IT_0015);
    const ccomplex_t IT_0017 = 0.666666666666667*IT_0016;
    const ccomplex_t IT_0018 = pow(m_t, 2);
    const ccomplex_t IT_0019 = pow(m_N_1, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*U_sb_01*conj(U_st_01);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0027 = IT_0005*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = m_t*N_u1*e_em*IT_0012*U_st_01;
    const ccomplex_t IT_0030 = IT_0011*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + (-0.75)*IT_0031);
    const ccomplex_t IT_0033 = 0.666666666666667*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = IT_0025*IT_0036;
    const ccomplex_t IT_0038 = pow(m_sb_2, 2);
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0018 + IT_0038 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0040 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0043 = IT_0001*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0041 + (-0.5)*IT_0044);
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0001*IT_0010*IT_0012*U_sb_01;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = IT_0039*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0018 + IT_0038 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0053 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0054 = IT_0001*IT_0053;
    const ccomplex_t IT_0055 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0056 = IT_0001*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0054 + (-0.5)*IT_0057);
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0001*IT_0010*IT_0012*U_sb_01;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = IT_0052*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = IT_0023 + IT_0037 + 2*IT_0051 + 2*IT_0064;
    const ccomplex_t IT_0066 = pow(m_W, -2);
    const ccomplex_t IT_0067 = pow(s_14, 2);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = -IT_0019;
    const ccomplex_t IT_0070 = IT_0068 + IT_0069;
    const ccomplex_t IT_0071 = s_13*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0073 = cos(beta);
    const ccomplex_t IT_0074 = cpow(IT_0073, -1);
    const ccomplex_t IT_0075 = IT_0001*IT_0074;
    const ccomplex_t IT_0076 = m_b*conj(U_d1)*V_tb*e_em*IT_0012*U_sb_11;
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0080 = IT_0001*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0078 + (-2)*IT_0080);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0046*IT_0082;
    const ccomplex_t IT_0084 = IT_0039*IT_0072*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0086 = m_b*conj(U_d2)*V_tb*e_em*IT_0012*U_sb_11;
    const ccomplex_t IT_0087 = IT_0075*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0090 = IT_0001*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0088 + (-2)*IT_0090);
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = IT_0059*IT_0092;
    const ccomplex_t IT_0094 = IT_0052*IT_0085*IT_0093;
    const ccomplex_t IT_0095 = pow(m_b, 2);
    const ccomplex_t IT_0096 = cpow((-2)*s_12 + -IT_0019 + -IT_0038 + IT_0095 
      + -reg_prop, -1);
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0098 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0099 = IT_0005*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = m_b*N_d1*e_em*IT_0012*U_sb_01;
    const ccomplex_t IT_0102 = IT_0075*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + 1.5*IT_0103);
    const ccomplex_t IT_0105 = (-0.333333333333333)*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = IT_0105*IT_0107;
    const ccomplex_t IT_0109 = IT_0096*IT_0097*IT_0108;
    const ccomplex_t IT_0110 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0111 = IT_0001*IT_0110;
    const ccomplex_t IT_0112 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0113 = IT_0001*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + 0.5*IT_0114);
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = IT_0092*IT_0116;
    const ccomplex_t IT_0118 = IT_0052*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = m_N_1*IT_0119;
    const ccomplex_t IT_0121 = -IT_0084 + -IT_0094 + -IT_0109 + IT_0120;
    const ccomplex_t IT_0122 = m_t*IT_0070;
    const ccomplex_t IT_0123 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0124 = IT_0005*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0127 = IT_0001*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = m_t*conj(N_u1)*e_em*IT_0012*U_st_11;
    const ccomplex_t IT_0130 = IT_0011*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0125 + 3*IT_0128 + 3
      *IT_0131);
    const ccomplex_t IT_0133 = 0.166666666666667*IT_0132;
    const ccomplex_t IT_0134 = IT_0034*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0025*IT_0135;
    const ccomplex_t IT_0137 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0138 = IT_0005*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0141 = IT_0001*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = m_t*conj(N_u1)*e_em*IT_0012*U_st_10;
    const ccomplex_t IT_0144 = IT_0011*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*(IT_0139 + 3*IT_0142 + 3
      *IT_0145);
    const ccomplex_t IT_0147 = 0.166666666666667*IT_0146;
    const ccomplex_t IT_0148 = IT_0020*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0003*IT_0149;
    const ccomplex_t IT_0151 = (-2)*IT_0119 + -IT_0136 + -IT_0150;
    const ccomplex_t IT_0152 = m_t*m_N_1;
    const ccomplex_t IT_0153 = IT_0070*IT_0152;
    const ccomplex_t IT_0154 = s_14*s_34*IT_0066;
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = s_13 + IT_0155;
    const ccomplex_t IT_0157 = s_13*IT_0156;
    const ccomplex_t IT_0158 = -IT_0023 + -IT_0037;
    const ccomplex_t IT_0159 = (-6)*IT_0065;
    const ccomplex_t IT_0160 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0161 = IT_0152*IT_0156;
    const ccomplex_t IT_0162 = pow(s_34, 2);
    const ccomplex_t IT_0163 = IT_0066*IT_0162;
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = IT_0018 + IT_0164;
    const ccomplex_t IT_0166 = s_13*IT_0165;
    const ccomplex_t IT_0167 = m_t*IT_0156;
    const ccomplex_t IT_0168 = (-6)*IT_0158;
    const ccomplex_t IT_0169 = m_N_1*IT_0156;
    const ccomplex_t IT_0170 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0171 = IT_0001*IT_0170;
    const ccomplex_t IT_0172 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0173 = IT_0001*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0171 + 0.5*IT_0174);
    const ccomplex_t IT_0176 = -IT_0175;
    const ccomplex_t IT_0177 = IT_0082*IT_0176;
    const ccomplex_t IT_0178 = IT_0039*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = pow(m_W, 2);
    const ccomplex_t IT_0182 = s_34*IT_0181;
    const ccomplex_t IT_0183 = IT_0066*IT_0182;
    const ccomplex_t IT_0184 = (-0.25)*IT_0183;
    const ccomplex_t IT_0185 = s_34 + IT_0184;
    const ccomplex_t IT_0186 = m_N_1*IT_0185;
    const ccomplex_t IT_0187 = (-4)*IT_0186;
    const ccomplex_t IT_0188 = (-6)*IT_0121;
    const ccomplex_t IT_0189 = m_b*conj(N_d1)*e_em*IT_0012*U_sb_11;
    const ccomplex_t IT_0190 = IT_0075*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0193 = IT_0001*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0196 = IT_0005*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*(IT_0191 + -IT_0194 +
       0.333333333333333*IT_0197);
    const ccomplex_t IT_0199 = 3*IT_0198;
    const ccomplex_t IT_0200 = 0.166666666666667*IT_0199;
    const ccomplex_t IT_0201 = IT_0107*IT_0200;
    const ccomplex_t IT_0202 = IT_0096*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = IT_0136 + IT_0150 + 2*IT_0203;
    const ccomplex_t IT_0205 = (-6)*IT_0151;
    const ccomplex_t IT_0206 = IT_0152*IT_0165;
    const ccomplex_t IT_0207 = (-6)*IT_0204;
    const ccomplex_t IT_0208 = m_N_1*IT_0165;
    const ccomplex_t IT_0209 = m_N_1*IT_0051;
    const ccomplex_t IT_0210 = m_N_1*IT_0064;
    const ccomplex_t IT_0211 = m_t*IT_0203;
    const ccomplex_t IT_0212 = -IT_0209 + -IT_0210 + -IT_0211;
    const ccomplex_t IT_0213 = IT_0048*IT_0176;
    const ccomplex_t IT_0214 = IT_0039*IT_0072*IT_0213;
    const ccomplex_t IT_0215 = IT_0061*IT_0116;
    const ccomplex_t IT_0216 = IT_0052*IT_0085*IT_0215;
    const ccomplex_t IT_0217 = IT_0214 + IT_0216;
    const ccomplex_t IT_0218 = m_N_1*IT_0179;
    const ccomplex_t IT_0219 = 2*IT_0179;
    const ccomplex_t IT_0220 = (-6)*IT_0219;
    const ccomplex_t IT_0221 = 0.166666666666667*IT_0220;
    const ccomplex_t IT_0222 = IT_0151 + IT_0221;
    const ccomplex_t IT_0223 = s_13*IT_0181;
    const ccomplex_t IT_0224 = s_14*s_34;
    const ccomplex_t IT_0225 = (-2)*IT_0224;
    const ccomplex_t IT_0226 = IT_0223 + IT_0225;
    const ccomplex_t IT_0227 = IT_0066*IT_0226;
    const ccomplex_t IT_0228 = (-6)*IT_0227;
    const ccomplex_t IT_0229 = 12*s_13;
    const ccomplex_t IT_0230 = IT_0228 + IT_0229;
    const ccomplex_t IT_0231 = IT_0212 + IT_0217;
    const ccomplex_t IT_0232 = -IT_0051 + -IT_0064;
    const ccomplex_t IT_0233 = 6*IT_0232;
    const ccomplex_t IT_0234 = s_14*IT_0181;
    const ccomplex_t IT_0235 = IT_0066*IT_0234;
    const ccomplex_t IT_0236 = (-0.25)*IT_0235;
    const ccomplex_t IT_0237 = s_14 + IT_0236;
    const ccomplex_t IT_0238 = m_t*IT_0237;
    const ccomplex_t IT_0239 = (-4)*IT_0238;
    const ccomplex_t IT_0240 = 6*IT_0203;
    const ccomplex_t IT_0241 = 0.5*IT_0235;
    const ccomplex_t IT_0242 = s_14 + IT_0241;
    const ccomplex_t IT_0243 = m_t*IT_0242;
    const ccomplex_t IT_0244 = 2*IT_0243;
    const ccomplex_t IT_0245 = 6*IT_0119 + (-6)*IT_0180;
    const ccomplex_t IT_0246 = IT_0159*IT_0169 + IT_0167*IT_0207 + IT_0168
      *IT_0208 + IT_0152*((-18)*IT_0121 + (-18)*IT_0218) + 6*IT_0122*IT_0222 +
       IT_0230*IT_0231 + IT_0187*IT_0233 + IT_0239*IT_0240 + IT_0244*IT_0245;
    const ccomplex_t IT_0247 = 0.166666666666667*IT_0205;
    const ccomplex_t IT_0248 = IT_0219 + IT_0247;
    const ccomplex_t IT_0249 = 0.5*IT_0183;
    const ccomplex_t IT_0250 = s_34 + IT_0249;
    const ccomplex_t IT_0251 = m_N_1*IT_0250;
    const ccomplex_t IT_0252 = 2*IT_0251;
    const ccomplex_t IT_0253 = 6*IT_0065*IT_0122 + IT_0167*IT_0168 + IT_0207
      *IT_0208 + IT_0152*((-18)*IT_0212 + (-18)*IT_0217) + (IT_0121 + IT_0218)
      *IT_0230 + IT_0233*IT_0244 + IT_0187*IT_0245 + 6*IT_0169*IT_0248 + IT_0240
      *IT_0252;
    const ccomplex_t IT_0254 = IT_0152*IT_0181;
    const ccomplex_t IT_0255 = 3*IT_0254;
    const ccomplex_t IT_0256 = pow(m_W, 4);
    const ccomplex_t IT_0257 = s_13*IT_0256;
    const ccomplex_t IT_0258 = IT_0066*IT_0257;
    const ccomplex_t IT_0259 = 24*IT_0224;
    const ccomplex_t IT_0260 = (-12)*IT_0223;
    const ccomplex_t IT_0261 = 6*IT_0258 + IT_0259 + IT_0260;
    const ccomplex_t IT_0262 = (-3)*IT_0254;
    const ccomplex_t IT_0263 = -conj(IT_0119);
    const ccomplex_t IT_0264 = conj(IT_0180) + IT_0263;
    const ccomplex_t IT_0265 = -IT_0264;
    const ccomplex_t IT_0266 = s_13*s_14;
    const ccomplex_t IT_0267 = s_34*IT_0019;
    const ccomplex_t IT_0268 = -IT_0267;
    const ccomplex_t IT_0269 = IT_0266 + IT_0268;
    const ccomplex_t IT_0270 = (-6)*conj(IT_0219);
    const ccomplex_t IT_0271 = 0.166666666666667*IT_0270;
    const ccomplex_t IT_0272 = conj(IT_0151) + IT_0271;
    const ccomplex_t IT_0273 = IT_0180 + IT_0203;
    const ccomplex_t IT_0274 = conj(IT_0180) + conj(IT_0203);
    const ccomplex_t IT_0275 = (-6)*conj(IT_0151);
    const ccomplex_t IT_0276 = 0.166666666666667*IT_0273;
    const ccomplex_t IT_0277 = s_14*IT_0018;
    const ccomplex_t IT_0278 = s_13*s_34;
    const ccomplex_t IT_0279 = -IT_0278;
    const ccomplex_t IT_0280 = IT_0277 + IT_0279;
    const ccomplex_t IT_0281 = (-6)*conj(IT_0204);
    const ccomplex_t IT_0282 = (-6)*IT_0258;
    const ccomplex_t IT_0283 = IT_0259 + IT_0282;
    const ccomplex_t IT_0284 = 0.166666666666667*conj(IT_0180);
    const ccomplex_t IT_0285 = (-24)*IT_0224;
    const ccomplex_t IT_0286 = 6*IT_0258 + IT_0285;
    const ccomplex_t IT_0287 = 0.166666666666667*conj(IT_0119);
    const ccomplex_t IT_0288 = 0.166666666666667*conj(IT_0203);
    const ccomplex_t IT_0289 = 0.166666666666667*IT_0153;
    const ccomplex_t IT_0290 = 0.166666666666667*IT_0275;
    const ccomplex_t IT_0291 = conj(IT_0219) + IT_0290;
    const ccomplex_t IT_0292 = (-6)*conj(IT_0158);
    const ccomplex_t IT_0293 = 0.166666666666667*IT_0161;
    const ccomplex_t IT_0294 = (-6)*conj(IT_0121);
    const ccomplex_t IT_0295 = 0.166666666666667*IT_0169;
    const ccomplex_t IT_0296 = conj(IT_0065)*IT_0122;
    const ccomplex_t IT_0297 = 0.166666666666667*IT_0167;
    const ccomplex_t IT_0298 = (-18)*conj(IT_0121);
    const ccomplex_t IT_0299 = 0.166666666666667*IT_0208;
    const ccomplex_t IT_0300 = 12*IT_0223;
    const ccomplex_t IT_0301 = IT_0282 + IT_0285 + IT_0300;
    const ccomplex_t IT_0302 = conj(IT_0232)*IT_0262;
    const ccomplex_t IT_0303 = 6*IT_0065*(conj(IT_0065)*IT_0071 + conj(IT_0121
      )*IT_0122 + conj(IT_0151)*IT_0153) + IT_0157*(conj(IT_0158)*IT_0159 +
       IT_0158*IT_0160) + (conj(IT_0151)*IT_0161 + conj(IT_0158)*IT_0166 + conj
      (IT_0121)*IT_0167)*IT_0168 + (conj(IT_0151)*IT_0169 + conj(IT_0180)
      *IT_0187)*IT_0188 + conj(IT_0204)*(IT_0159*IT_0161 + IT_0157*IT_0205 +
       IT_0168*IT_0206 + IT_0166*IT_0207 + IT_0188*IT_0208) + (conj(IT_0212) +
       conj(IT_0217))*IT_0246 + conj(IT_0218)*IT_0253 + IT_0233*(conj(IT_0121)
      *IT_0244 + conj(IT_0203)*IT_0255 + 0.166666666666667*conj(IT_0232)*IT_0261
       + IT_0262*IT_0265) + 6*IT_0269*(conj(IT_0119)*IT_0222 + conj(IT_0065)
      *IT_0232 + IT_0065*conj(IT_0232) + IT_0119*IT_0272 + conj(IT_0219)*IT_0273
       + IT_0248*IT_0274 + IT_0275*IT_0276) + 6*IT_0280*(conj(IT_0119)*IT_0204 +
       IT_0119*conj(IT_0204) + conj(IT_0158)*IT_0232 + IT_0158*conj(IT_0232) +
       0.166666666666667*IT_0207*IT_0274 + IT_0276*IT_0281) + IT_0240*(conj
      (IT_0121)*IT_0252 + conj(IT_0232)*IT_0255 + IT_0283*IT_0284 + IT_0286
      *IT_0287 + IT_0261*IT_0288) + 6*conj(IT_0219)*(IT_0158*IT_0161 +
       0.166666666666667*IT_0071*IT_0205 + IT_0159*IT_0289) + 6*IT_0219*(conj
      (IT_0158)*IT_0161 + conj(IT_0121)*IT_0169 + IT_0157*conj(IT_0204) +
       IT_0160*IT_0289 + IT_0071*IT_0291) + 6*IT_0151*(IT_0071*conj(IT_0151) +
       conj(IT_0065)*IT_0153 + IT_0292*IT_0293 + IT_0294*IT_0295) + 6*IT_0121*
      (conj(IT_0119)*IT_0187 + IT_0169*conj(IT_0219) + 0.166666666666667*conj
      (IT_0121)*IT_0230 + conj(IT_0232)*IT_0244 + conj(IT_0203)*IT_0252 +
       IT_0296 + IT_0292*IT_0297) + 6*IT_0218*(0.166666666666667*conj(IT_0121)
      *IT_0230 + conj(IT_0232)*IT_0244 + conj(IT_0203)*IT_0252 + IT_0187*IT_0265
       + IT_0169*IT_0291 + IT_0296 + IT_0292*IT_0297 + IT_0281*IT_0299) + 6
      *IT_0231*(IT_0187*conj(IT_0232) + conj(IT_0203)*IT_0239 + IT_0244*IT_0265 
      + IT_0122*IT_0272 + IT_0160*IT_0295 + IT_0281*IT_0297 + 0.166666666666667
      *IT_0152*IT_0298 + IT_0292*IT_0299) + 6*IT_0204*(IT_0157*IT_0291 +
       0.166666666666667*IT_0206*IT_0292 + IT_0160*IT_0293 + IT_0294*IT_0299) +
       6*IT_0119*(conj(IT_0121)*IT_0187 + IT_0261*IT_0287 + IT_0286*IT_0288 +
       IT_0284*IT_0301 + IT_0302) + IT_0180*(conj(IT_0180)*IT_0261 + conj
      (IT_0203)*IT_0283 + IT_0187*IT_0294 + conj(IT_0119)*IT_0301 + (-6)*IT_0302);
    return create_ccomplex_return(IT_0303);
}

