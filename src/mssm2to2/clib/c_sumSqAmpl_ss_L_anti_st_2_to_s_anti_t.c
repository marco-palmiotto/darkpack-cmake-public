#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_st_2_to_s_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_st_2_to_s_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
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
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B2)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W2)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + (-3)*IT_0009);
    const ccomplex_t IT_0011 = 0.166666666666667*IT_0010;
    const ccomplex_t IT_0012 = N_B2*e_em*conj(U_st_01);
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = N_W2*e_em*conj(U_st_01);
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = sin(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0006*IT_0019;
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = m_t*N_u2*e_em*IT_0021*conj(U_st_11);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0014 + 3*IT_0017 + 3
      *IT_0024);
    const ccomplex_t IT_0026 = 0.166666666666667*IT_0025;
    const ccomplex_t IT_0027 = IT_0011*IT_0026;
    const ccomplex_t IT_0028 = pow(m_s, 2);
    const ccomplex_t IT_0029 = pow(m_ss_L, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0028 + IT_0029 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0027*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = conj(N_B3)*e_em;
    const ccomplex_t IT_0034 = IT_0001*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = conj(N_W3)*e_em;
    const ccomplex_t IT_0037 = IT_0006*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + (-3)*IT_0038);
    const ccomplex_t IT_0040 = 0.166666666666667*IT_0039;
    const ccomplex_t IT_0041 = N_B3*e_em*conj(U_st_01);
    const ccomplex_t IT_0042 = IT_0001*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = N_W3*e_em*conj(U_st_01);
    const ccomplex_t IT_0045 = IT_0006*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = m_t*N_u3*e_em*IT_0021*conj(U_st_11);
    const ccomplex_t IT_0048 = IT_0020*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0043 + 3*IT_0046 + 3
      *IT_0049);
    const ccomplex_t IT_0051 = 0.166666666666667*IT_0050;
    const ccomplex_t IT_0052 = IT_0040*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_13 + IT_0028 + IT_0029 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = conj(N_B4)*e_em;
    const ccomplex_t IT_0057 = IT_0001*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W4)*e_em;
    const ccomplex_t IT_0060 = IT_0006*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + (-3)*IT_0061);
    const ccomplex_t IT_0063 = 0.166666666666667*IT_0062;
    const ccomplex_t IT_0064 = N_B4*e_em*conj(U_st_01);
    const ccomplex_t IT_0065 = IT_0001*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = N_W4*e_em*conj(U_st_01);
    const ccomplex_t IT_0068 = IT_0006*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = m_t*N_u4*e_em*IT_0021*conj(U_st_11);
    const ccomplex_t IT_0071 = IT_0020*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0066 + 3*IT_0069 + 3
      *IT_0072);
    const ccomplex_t IT_0074 = 0.166666666666667*IT_0073;
    const ccomplex_t IT_0075 = IT_0063*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_13 + IT_0028 + IT_0029 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_01);
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = pow(m_sG, 2);
    const ccomplex_t IT_0083 = cpow((-2)*s_13 + IT_0028 + IT_0029 + -IT_0082 +
       reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_13 + IT_0028 + IT_0029 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0087 = N_B1*e_em*conj(U_st_01);
    const ccomplex_t IT_0088 = IT_0001*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = N_W1*e_em*conj(U_st_01);
    const ccomplex_t IT_0091 = IT_0006*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = m_t*N_u1*e_em*IT_0021*conj(U_st_11);
    const ccomplex_t IT_0094 = IT_0020*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0089 + 3*IT_0092 + 3
      *IT_0095);
    const ccomplex_t IT_0097 = 0.166666666666667*IT_0096;
    const ccomplex_t IT_0098 = conj(N_B1)*e_em;
    const ccomplex_t IT_0099 = IT_0001*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W1)*e_em;
    const ccomplex_t IT_0102 = IT_0006*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + (-3)*IT_0103);
    const ccomplex_t IT_0105 = 0.166666666666667*IT_0104;
    const ccomplex_t IT_0106 = IT_0097*IT_0105;
    const ccomplex_t IT_0107 = IT_0086*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = -IT_0032 + -IT_0055 + -IT_0078 +
       0.166666666666667*IT_0085 + -IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0006*conj(U_st_01);
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0006;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = pow(m_st_2, 2);
    const ccomplex_t IT_0115 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0029 + -IT_0114 + -reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0113*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = IT_0111*IT_0117;
    const ccomplex_t IT_0119 = (-0.5)*IT_0085 + -IT_0118;
    const ccomplex_t IT_0120 = cos(beta);
    const ccomplex_t IT_0121 = cpow(IT_0120, -1);
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0006*IT_0021*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = conj(N_B2)*e_em*conj(U_st_11);
    const ccomplex_t IT_0125 = IT_0001*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = m_t*conj(N_u2)*e_em*IT_0021*conj(U_st_01);
    const ccomplex_t IT_0128 = IT_0020*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + (-0.75)*IT_0129);
    const ccomplex_t IT_0131 = 0.666666666666667*IT_0130;
    const ccomplex_t IT_0132 = IT_0123*IT_0131;
    const ccomplex_t IT_0133 = IT_0030*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0006*IT_0021*IT_0121;
    const ccomplex_t IT_0136 = (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = conj(N_B3)*e_em*conj(U_st_11);
    const ccomplex_t IT_0138 = IT_0001*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = m_t*conj(N_u3)*e_em*IT_0021*conj(U_st_01);
    const ccomplex_t IT_0141 = IT_0020*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + (-0.75)*IT_0142);
    const ccomplex_t IT_0144 = 0.666666666666667*IT_0143;
    const ccomplex_t IT_0145 = IT_0136*IT_0144;
    const ccomplex_t IT_0146 = IT_0053*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0006*IT_0021*IT_0121;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = conj(N_B4)*e_em*conj(U_st_11);
    const ccomplex_t IT_0151 = IT_0001*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = m_t*conj(N_u4)*e_em*IT_0021*conj(U_st_01);
    const ccomplex_t IT_0154 = IT_0020*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0152 + (-0.75)*IT_0155);
    const ccomplex_t IT_0157 = 0.666666666666667*IT_0156;
    const ccomplex_t IT_0158 = IT_0149*IT_0157;
    const ccomplex_t IT_0159 = IT_0076*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = conj(N_B1)*e_em*conj(U_st_11);
    const ccomplex_t IT_0162 = IT_0001*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = m_t*conj(N_u1)*e_em*IT_0021*conj(U_st_01);
    const ccomplex_t IT_0165 = IT_0020*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0163 + (-0.75)*IT_0166);
    const ccomplex_t IT_0168 = 0.666666666666667*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0006*IT_0021*IT_0121;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = IT_0168*IT_0170;
    const ccomplex_t IT_0172 = IT_0086*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = -IT_0134 + -IT_0147 + -IT_0160 + -IT_0173;
    const ccomplex_t IT_0175 = m_s*m_t*IT_0029;
    const ccomplex_t IT_0176 = s_34*IT_0029;
    const ccomplex_t IT_0177 = (-18)*IT_0176;
    const ccomplex_t IT_0178 = s_13*s_14;
    const ccomplex_t IT_0179 = 36*IT_0178;
    const ccomplex_t IT_0180 = IT_0177 + IT_0179;
    const ccomplex_t IT_0181 = m_s*s_14;
    const ccomplex_t IT_0182 = IT_0026*IT_0123;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0184 = IT_0030*IT_0182*IT_0183;
    const ccomplex_t IT_0185 = IT_0051*IT_0136;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0187 = IT_0053*IT_0185*IT_0186;
    const ccomplex_t IT_0188 = IT_0074*IT_0149;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0190 = IT_0076*IT_0188*IT_0189;
    const ccomplex_t IT_0191 = m_s*IT_0032;
    const ccomplex_t IT_0192 = m_s*IT_0078;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0194 = IT_0097*IT_0170;
    const ccomplex_t IT_0195 = IT_0086*IT_0193*IT_0194;
    const ccomplex_t IT_0196 = m_s*IT_0085;
    const ccomplex_t IT_0197 = m_s*IT_0108;
    const ccomplex_t IT_0198 = m_s*IT_0055;
    const ccomplex_t IT_0199 = -IT_0184 + -IT_0187 + -IT_0190 + IT_0191 +
       IT_0192 + -IT_0195 + (-0.166666666666667)*IT_0196 + IT_0197 + IT_0198;
    const ccomplex_t IT_0200 = A_t*IT_0120*conj(U_st_11) + (-1.4142135623731)
      *e_em*IT_0006*(m_W*IT_0018*IT_0120*conj(U_st_01) + (-0.5)*IT_0021*(IT_0018
      *IT_0028*IT_0121*conj(U_st_01) + m_t*(m_t*IT_0019*IT_0120*conj(U_st_01) +
       mu_h*conj(U_st_11))));
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*V_ts*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0006*IT_0018*IT_0021*IT_0121;
    const ccomplex_t IT_0203 = 0.5*IT_0202;
    const ccomplex_t IT_0204 = IT_0201*IT_0203;
    const ccomplex_t IT_0205 = cpow((-2)*s_12 + -IT_0029 + -IT_0114 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0206 = IT_0204*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*IT_0206;
    const ccomplex_t IT_0208 = (-0.5)*IT_0110;
    const ccomplex_t IT_0209 = m_s*IT_0208;
    const ccomplex_t IT_0210 = m_s*IT_0111;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0209 + IT_0211;
    const ccomplex_t IT_0213 = IT_0117*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0006*IT_0021;
    const ccomplex_t IT_0215 = (-0.5)*IT_0214;
    const ccomplex_t IT_0216 = pow(m_t, 2);
    const ccomplex_t IT_0217 = cpow(IT_0120, 2);
    const ccomplex_t IT_0218 = cpow(IT_0018, 2);
    const ccomplex_t IT_0219 = A_t*IT_0018*conj(U_st_11) + (-0.707106781186548
      )*e_em*IT_0006*((IT_0021*(IT_0028 + -IT_0216) + -m_W*(IT_0217 + -IT_0218))
      *conj(U_st_01) + m_t*mu_h*IT_0019*IT_0021*IT_0120*conj(U_st_11));
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*V_ts*IT_0219;
    const ccomplex_t IT_0221 = IT_0215*IT_0220;
    const ccomplex_t IT_0222 = pow(m_W, 2);
    const ccomplex_t IT_0223 = cpow((-2)*s_12 + -IT_0029 + -IT_0114 + IT_0222 
      + -reg_prop, -1);
    const ccomplex_t IT_0224 = IT_0221*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = 0.5*IT_0196 + -IT_0207 + -IT_0213 + -IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_11);
    const ccomplex_t IT_0228 = IT_0079*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0230 = IT_0083*IT_0228*IT_0229;
    const ccomplex_t IT_0231 = IT_0011*IT_0131;
    const ccomplex_t IT_0232 = IT_0030*IT_0183*IT_0231;
    const ccomplex_t IT_0233 = IT_0040*IT_0144;
    const ccomplex_t IT_0234 = IT_0053*IT_0186*IT_0233;
    const ccomplex_t IT_0235 = IT_0063*IT_0157;
    const ccomplex_t IT_0236 = IT_0076*IT_0189*IT_0235;
    const ccomplex_t IT_0237 = m_s*IT_0134;
    const ccomplex_t IT_0238 = m_s*IT_0147;
    const ccomplex_t IT_0239 = m_s*IT_0160;
    const ccomplex_t IT_0240 = IT_0105*IT_0168;
    const ccomplex_t IT_0241 = IT_0086*IT_0193*IT_0240;
    const ccomplex_t IT_0242 = m_s*IT_0173;
    const ccomplex_t IT_0243 = 0.166666666666667*IT_0230 + -IT_0232 + -IT_0234
       + -IT_0236 + IT_0237 + IT_0238 + IT_0239 + -IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0006*IT_0021;
    const ccomplex_t IT_0245 = 0.5*IT_0244;
    const ccomplex_t IT_0246 = IT_0220*IT_0245;
    const ccomplex_t IT_0247 = IT_0223*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = m_t*IT_0208;
    const ccomplex_t IT_0250 = m_t*IT_0111;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = IT_0249 + IT_0251;
    const ccomplex_t IT_0253 = IT_0117*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0006*IT_0019*IT_0021*IT_0120;
    const ccomplex_t IT_0255 = 0.5*IT_0254;
    const ccomplex_t IT_0256 = IT_0201*IT_0255;
    const ccomplex_t IT_0257 = IT_0205*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = (-0.5)*IT_0230 + -IT_0248 + IT_0253 + -IT_0258;
    const ccomplex_t IT_0260 = conj(IT_0199) + 3*conj(IT_0226);
    const ccomplex_t IT_0261 = 6*conj(IT_0119)*(IT_0199 + 3*IT_0226) + conj
      (IT_0109)*(18*IT_0199 + 6*IT_0226) + IT_0109*(18*conj(IT_0199) + 6*conj
      (IT_0226)) + conj(IT_0174)*(18*IT_0243 + 6*IT_0259) + IT_0174*(18*conj
      (IT_0243) + 6*conj(IT_0259)) + 6*IT_0119*IT_0260;
    const ccomplex_t IT_0262 = 6*s_34;
    const ccomplex_t IT_0263 = 18*s_34;
    const ccomplex_t IT_0264 = m_s*m_t;
    const ccomplex_t IT_0265 = (-6)*conj(IT_0259);
    const ccomplex_t IT_0266 = (-18)*conj(IT_0243);
    const ccomplex_t IT_0267 = IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = (-18)*conj(IT_0199);
    const ccomplex_t IT_0269 = (-6)*conj(IT_0226);
    const ccomplex_t IT_0270 = IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = (-6)*IT_0226;
    const ccomplex_t IT_0272 = m_t*s_13;
    const ccomplex_t IT_0273 = (-6)*conj(IT_0119)*(IT_0243 + 3*IT_0259) + conj
      (IT_0109)*((-18)*IT_0243 + (-6)*IT_0259) + (-6)*IT_0119*(conj(IT_0243) + 3
      *conj(IT_0259)) + IT_0109*IT_0267 + IT_0174*IT_0270 + conj(IT_0174)*((-18)
      *IT_0199 + IT_0271);
    const ccomplex_t IT_0274 = 12*IT_0178;
    const ccomplex_t IT_0275 = (-6)*IT_0176;
    const ccomplex_t IT_0276 = (((-18)*conj(IT_0109) + (-6)*conj(IT_0119))
      *IT_0174 + ((-18)*IT_0109 + (-6)*IT_0119)*conj(IT_0174))*IT_0175 + 
      (IT_0109*conj(IT_0109) + IT_0119*conj(IT_0119) + IT_0174*conj(IT_0174))
      *IT_0180 + IT_0181*IT_0261 + (conj(IT_0199)*IT_0226 + IT_0199*conj(IT_0226
      ) + conj(IT_0243)*IT_0259 + IT_0243*conj(IT_0259))*IT_0262 + (IT_0199*conj
      (IT_0199) + IT_0226*conj(IT_0226) + IT_0243*conj(IT_0243) + IT_0259*conj
      (IT_0259))*IT_0263 + IT_0264*((-18)*IT_0226*conj(IT_0259) + (-6)*IT_0259
      *IT_0260 + IT_0199*IT_0267 + IT_0243*IT_0270 + conj(IT_0243)*IT_0271) +
       IT_0272*IT_0273 + (conj(IT_0109)*IT_0119 + IT_0109*conj(IT_0119))*
      (IT_0274 + IT_0275);
    return create_ccomplex_return(IT_0276);
}

