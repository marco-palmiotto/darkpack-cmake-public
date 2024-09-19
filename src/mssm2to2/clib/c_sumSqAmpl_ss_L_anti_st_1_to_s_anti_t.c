#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_st_1_to_s_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_st_1_to_s_anti_t(
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
    const creal_t m_st_1 = param->m_st_1;
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
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_ss_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = m_t*N_u1*e_em*IT_0016*conj(U_st_10);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0007 + 3*IT_0012 + 3
      *IT_0019);
    const ccomplex_t IT_0021 = 0.166666666666667*IT_0020;
    const ccomplex_t IT_0022 = conj(N_B1)*e_em;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W1)*e_em;
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + (-3)*IT_0027);
    const ccomplex_t IT_0029 = 0.166666666666667*IT_0028;
    const ccomplex_t IT_0030 = IT_0021*IT_0029;
    const ccomplex_t IT_0031 = IT_0002*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0034 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0038 = IT_0009*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = m_t*N_u2*e_em*IT_0016*conj(U_st_10);
    const ccomplex_t IT_0041 = IT_0015*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0036 + 3*IT_0039 + 3
      *IT_0042);
    const ccomplex_t IT_0044 = 0.166666666666667*IT_0043;
    const ccomplex_t IT_0045 = conj(N_B2)*e_em;
    const ccomplex_t IT_0046 = IT_0004*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = conj(N_W2)*e_em;
    const ccomplex_t IT_0049 = IT_0009*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0047 + (-3)*IT_0050);
    const ccomplex_t IT_0052 = 0.166666666666667*IT_0051;
    const ccomplex_t IT_0053 = IT_0044*IT_0052;
    const ccomplex_t IT_0054 = IT_0033*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = conj(N_B3)*e_em;
    const ccomplex_t IT_0057 = IT_0004*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W3)*e_em;
    const ccomplex_t IT_0060 = IT_0009*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + (-3)*IT_0061);
    const ccomplex_t IT_0063 = 0.166666666666667*IT_0062;
    const ccomplex_t IT_0064 = N_B3*e_em*conj(U_st_00);
    const ccomplex_t IT_0065 = IT_0004*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = N_W3*e_em*conj(U_st_00);
    const ccomplex_t IT_0068 = IT_0009*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = m_t*N_u3*e_em*IT_0016*conj(U_st_10);
    const ccomplex_t IT_0071 = IT_0015*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0066 + 3*IT_0069 + 3
      *IT_0072);
    const ccomplex_t IT_0074 = 0.166666666666667*IT_0073;
    const ccomplex_t IT_0075 = IT_0063*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0080 = conj(N_B4)*e_em;
    const ccomplex_t IT_0081 = IT_0004*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = conj(N_W4)*e_em;
    const ccomplex_t IT_0084 = IT_0009*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + (-3)*IT_0085);
    const ccomplex_t IT_0087 = 0.166666666666667*IT_0086;
    const ccomplex_t IT_0088 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0089 = IT_0009*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = m_t*N_u4*e_em*IT_0016*conj(U_st_10);
    const ccomplex_t IT_0092 = IT_0015*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0095 = IT_0004*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0090 + IT_0093 +
       0.333333333333333*IT_0096);
    const ccomplex_t IT_0098 = 3*IT_0097;
    const ccomplex_t IT_0099 = 0.166666666666667*IT_0098;
    const ccomplex_t IT_0100 = IT_0087*IT_0099;
    const ccomplex_t IT_0101 = IT_0079*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = -IT_0032 + -IT_0055 + -IT_0078 + -IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0009*conj(U_st_00);
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0009;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = pow(m_st_1, 2);
    const ccomplex_t IT_0109 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0001 + -IT_0108 + -reg_prop, -1);
    const ccomplex_t IT_0110 = IT_0107*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0105*IT_0111;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_00);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = pow(m_sG, 2);
    const ccomplex_t IT_0118 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0117 +
       reg_prop, -1);
    const ccomplex_t IT_0119 = IT_0116*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = 0.166666666666667*IT_0120;
    const ccomplex_t IT_0123 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0124 = IT_0004*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = m_t*conj(N_u1)*e_em*IT_0016*conj(U_st_00);
    const ccomplex_t IT_0127 = IT_0015*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0125 + (-0.75)*IT_0128);
    const ccomplex_t IT_0130 = 0.666666666666667*IT_0129;
    const ccomplex_t IT_0131 = cos(beta);
    const ccomplex_t IT_0132 = cpow(IT_0131, -1);
    const ccomplex_t IT_0133 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0009*IT_0016*IT_0132;
    const ccomplex_t IT_0134 = (-0.5)*IT_0133;
    const ccomplex_t IT_0135 = IT_0130*IT_0134;
    const ccomplex_t IT_0136 = IT_0002*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0009*IT_0016*IT_0132;
    const ccomplex_t IT_0139 = (-0.5)*IT_0138;
    const ccomplex_t IT_0140 = conj(N_B3)*e_em*conj(U_st_10);
    const ccomplex_t IT_0141 = IT_0004*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = m_t*conj(N_u3)*e_em*IT_0016*conj(U_st_00);
    const ccomplex_t IT_0144 = IT_0015*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*(IT_0142 + (-0.75)*IT_0145);
    const ccomplex_t IT_0147 = 0.666666666666667*IT_0146;
    const ccomplex_t IT_0148 = IT_0139*IT_0147;
    const ccomplex_t IT_0149 = IT_0076*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0009*IT_0016*IT_0132;
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0154 = IT_0004*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = m_t*conj(N_u4)*e_em*IT_0016*conj(U_st_00);
    const ccomplex_t IT_0157 = IT_0015*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0155 + (-0.75)*IT_0158);
    const ccomplex_t IT_0160 = 0.666666666666667*IT_0159;
    const ccomplex_t IT_0161 = IT_0152*IT_0160;
    const ccomplex_t IT_0162 = IT_0079*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0165 = IT_0004*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = m_t*conj(N_u2)*e_em*IT_0016*conj(U_st_00);
    const ccomplex_t IT_0168 = IT_0015*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*(IT_0166 + (-0.75)*IT_0169);
    const ccomplex_t IT_0171 = 0.666666666666667*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0009*IT_0016*IT_0132;
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = IT_0171*IT_0173;
    const ccomplex_t IT_0175 = IT_0033*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = -IT_0137 + -IT_0150 + -IT_0163 + -IT_0176;
    const ccomplex_t IT_0178 = m_s*m_t*IT_0001;
    const ccomplex_t IT_0179 = s_34*IT_0001;
    const ccomplex_t IT_0180 = (-6)*IT_0179;
    const ccomplex_t IT_0181 = s_13*s_14;
    const ccomplex_t IT_0182 = 12*IT_0181;
    const ccomplex_t IT_0183 = IT_0180 + IT_0182;
    const ccomplex_t IT_0184 = conj(IT_0103) + conj(IT_0122);
    const ccomplex_t IT_0185 = IT_0113 + IT_0121;
    const ccomplex_t IT_0186 = IT_0103 + IT_0122;
    const ccomplex_t IT_0187 = conj(IT_0113) + conj(IT_0121);
    const ccomplex_t IT_0188 = (-18)*IT_0179;
    const ccomplex_t IT_0189 = 36*IT_0181;
    const ccomplex_t IT_0190 = IT_0188 + IT_0189;
    const ccomplex_t IT_0191 = IT_0074*IT_0139;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0193 = IT_0076*IT_0191*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0195 = IT_0099*IT_0152;
    const ccomplex_t IT_0196 = IT_0079*IT_0194*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0198 = IT_0021*IT_0134;
    const ccomplex_t IT_0199 = IT_0002*IT_0197*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0201 = IT_0044*IT_0173;
    const ccomplex_t IT_0202 = IT_0033*IT_0200*IT_0201;
    const ccomplex_t IT_0203 = m_s*IT_0120;
    const ccomplex_t IT_0204 = m_s*IT_0032;
    const ccomplex_t IT_0205 = m_s*IT_0055;
    const ccomplex_t IT_0206 = m_s*IT_0078;
    const ccomplex_t IT_0207 = m_s*IT_0102;
    const ccomplex_t IT_0208 = -IT_0193 + -IT_0196 + -IT_0199 + -IT_0202 + (
      -0.166666666666667)*IT_0203 + IT_0204 + IT_0205 + IT_0206 + IT_0207;
    const ccomplex_t IT_0209 = A_t*IT_0131*conj(U_st_10) + (-1.4142135623731)
      *e_em*IT_0009*(m_W*IT_0013*IT_0131*conj(U_st_00) + (-0.5)*IT_0016*(IT_0000
      *IT_0013*IT_0132*conj(U_st_00) + m_t*(m_t*IT_0014*IT_0131*conj(U_st_00) +
       mu_h*conj(U_st_10))));
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*V_ts*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0009*IT_0013*IT_0016*IT_0132;
    const ccomplex_t IT_0212 = 0.5*IT_0211;
    const ccomplex_t IT_0213 = IT_0210*IT_0212;
    const ccomplex_t IT_0214 = cpow((-2)*s_12 + -IT_0001 + -IT_0108 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0215 = IT_0213*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = (-0.5)*IT_0104;
    const ccomplex_t IT_0218 = m_s*IT_0217;
    const ccomplex_t IT_0219 = m_s*IT_0105;
    const ccomplex_t IT_0220 = -IT_0219;
    const ccomplex_t IT_0221 = IT_0218 + IT_0220;
    const ccomplex_t IT_0222 = IT_0111*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0009*IT_0016;
    const ccomplex_t IT_0224 = (-0.5)*IT_0223;
    const ccomplex_t IT_0225 = pow(m_t, 2);
    const ccomplex_t IT_0226 = cpow(IT_0131, 2);
    const ccomplex_t IT_0227 = cpow(IT_0013, 2);
    const ccomplex_t IT_0228 = A_t*IT_0013*conj(U_st_10) + (-0.707106781186548
      )*e_em*IT_0009*((IT_0016*(IT_0000 + -IT_0225) + -m_W*(IT_0226 + -IT_0227))
      *conj(U_st_00) + m_t*mu_h*IT_0014*IT_0016*IT_0131*conj(U_st_10));
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*V_ts*IT_0228;
    const ccomplex_t IT_0230 = IT_0224*IT_0229;
    const ccomplex_t IT_0231 = pow(m_W, 2);
    const ccomplex_t IT_0232 = cpow((-2)*s_12 + -IT_0001 + -IT_0108 + IT_0231 
      + -reg_prop, -1);
    const ccomplex_t IT_0233 = IT_0230*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = 0.5*IT_0203 + -IT_0216 + -IT_0222 + -IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_10);
    const ccomplex_t IT_0237 = IT_0114*IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0239 = IT_0118*IT_0237*IT_0238;
    const ccomplex_t IT_0240 = IT_0063*IT_0147;
    const ccomplex_t IT_0241 = IT_0076*IT_0192*IT_0240;
    const ccomplex_t IT_0242 = IT_0087*IT_0160;
    const ccomplex_t IT_0243 = IT_0079*IT_0194*IT_0242;
    const ccomplex_t IT_0244 = m_s*IT_0150;
    const ccomplex_t IT_0245 = m_s*IT_0163;
    const ccomplex_t IT_0246 = IT_0029*IT_0130;
    const ccomplex_t IT_0247 = IT_0002*IT_0197*IT_0246;
    const ccomplex_t IT_0248 = IT_0052*IT_0171;
    const ccomplex_t IT_0249 = IT_0033*IT_0200*IT_0248;
    const ccomplex_t IT_0250 = m_s*IT_0176;
    const ccomplex_t IT_0251 = m_s*IT_0137;
    const ccomplex_t IT_0252 = 0.166666666666667*IT_0239 + -IT_0241 + -IT_0243
       + IT_0244 + IT_0245 + -IT_0247 + -IT_0249 + IT_0250 + IT_0251;
    const ccomplex_t IT_0253 = (-0.5)*IT_0239;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0009*IT_0014*IT_0016*IT_0131;
    const ccomplex_t IT_0255 = 0.5*IT_0254;
    const ccomplex_t IT_0256 = IT_0210*IT_0255;
    const ccomplex_t IT_0257 = IT_0214*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = m_t*IT_0217;
    const ccomplex_t IT_0260 = m_t*IT_0105;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0259 + IT_0261;
    const ccomplex_t IT_0263 = IT_0111*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0009*IT_0016;
    const ccomplex_t IT_0265 = 0.5*IT_0264;
    const ccomplex_t IT_0266 = IT_0229*IT_0265;
    const ccomplex_t IT_0267 = IT_0232*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*IT_0267;
    const ccomplex_t IT_0269 = -IT_0258 + IT_0263 + -IT_0268;
    const ccomplex_t IT_0270 = IT_0253 + IT_0269;
    const ccomplex_t IT_0271 = conj(IT_0253) + conj(IT_0269);
    const ccomplex_t IT_0272 = 6*s_34;
    const ccomplex_t IT_0273 = 18*s_34;
    const ccomplex_t IT_0274 = m_s*s_14;
    const ccomplex_t IT_0275 = 6*IT_0235;
    const ccomplex_t IT_0276 = 18*IT_0208;
    const ccomplex_t IT_0277 = 6*conj(IT_0235);
    const ccomplex_t IT_0278 = 18*conj(IT_0208);
    const ccomplex_t IT_0279 = 6*IT_0208;
    const ccomplex_t IT_0280 = 18*IT_0235;
    const ccomplex_t IT_0281 = 6*conj(IT_0208);
    const ccomplex_t IT_0282 = 18*conj(IT_0235);
    const ccomplex_t IT_0283 = 3*conj(IT_0177)*(IT_0252 + 1./3*IT_0253 + 1./3
      *IT_0269) + 3*IT_0177*(conj(IT_0252) + 1./3*conj(IT_0253) + 1./3*conj
      (IT_0269)) + 0.166666666666667*IT_0184*(IT_0275 + IT_0276) +
       0.166666666666667*IT_0186*(IT_0277 + IT_0278) + 0.166666666666667*IT_0187
      *(IT_0279 + IT_0280) + 0.166666666666667*IT_0185*(IT_0281 + IT_0282);
    const ccomplex_t IT_0284 = m_s*m_t;
    const ccomplex_t IT_0285 = (-18)*IT_0208;
    const ccomplex_t IT_0286 = (-6)*IT_0235;
    const ccomplex_t IT_0287 = IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = (-18)*conj(IT_0208);
    const ccomplex_t IT_0289 = (-6)*conj(IT_0235);
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = (-18)*IT_0235;
    const ccomplex_t IT_0292 = (-18)*conj(IT_0235);
    const ccomplex_t IT_0293 = (-6)*conj(IT_0208);
    const ccomplex_t IT_0294 = m_t*s_13;
    const ccomplex_t IT_0295 = (-18)*IT_0269;
    const ccomplex_t IT_0296 = (-6)*IT_0252;
    const ccomplex_t IT_0297 = (-18)*conj(IT_0253);
    const ccomplex_t IT_0298 = (-18)*conj(IT_0269);
    const ccomplex_t IT_0299 = (-6)*conj(IT_0252);
    const ccomplex_t IT_0300 = (-18)*IT_0252;
    const ccomplex_t IT_0301 = (-6)*IT_0269;
    const ccomplex_t IT_0302 = (-18)*conj(IT_0252);
    const ccomplex_t IT_0303 = (-6)*conj(IT_0253);
    const ccomplex_t IT_0304 = (-6)*conj(IT_0269);
    const ccomplex_t IT_0305 = (conj(IT_0103) + 3*conj(IT_0113) + 3*conj
      (IT_0121) + conj(IT_0122))*IT_0253 + (-0.166666666666667)*conj(IT_0177)
      *IT_0287 + (-0.166666666666667)*IT_0177*IT_0290 + (-0.166666666666667)
      *IT_0187*(IT_0295 + IT_0296) + (-0.166666666666667)*IT_0185*(IT_0297 +
       IT_0298 + IT_0299) + (-0.166666666666667)*IT_0184*(IT_0300 + IT_0301) + (
      -0.166666666666667)*IT_0186*(IT_0302 + IT_0303 + IT_0304);
    const ccomplex_t IT_0306 = (((-18)*conj(IT_0103) + (-6)*conj(IT_0113) + (
      -6)*conj(IT_0121) + (-18)*conj(IT_0122))*IT_0177 + ((-18)*IT_0103 + (-6)
      *IT_0113 + (-6)*IT_0121 + (-18)*IT_0122)*conj(IT_0177))*IT_0178 + IT_0183*
      (IT_0184*IT_0185 + IT_0186*IT_0187) + (IT_0177*conj(IT_0177) + (IT_0103 +
       IT_0122)*IT_0184 + IT_0185*IT_0187)*IT_0190 + (conj(IT_0208)*IT_0235 +
       IT_0208*conj(IT_0235) + conj(IT_0252)*IT_0270 + IT_0252*IT_0271)*IT_0272 
      + (IT_0208*conj(IT_0208) + IT_0235*conj(IT_0235) + IT_0252*conj(IT_0252) +
       IT_0270*IT_0271)*IT_0273 + 6*IT_0274*IT_0283 + IT_0284*(conj(IT_0252)
      *IT_0287 + IT_0252*IT_0290 + (-6)*IT_0271*(IT_0208 + (-0.166666666666667)
      *IT_0291) + IT_0270*(IT_0292 + IT_0293)) + (-6)*IT_0294*IT_0305;
    return create_ccomplex_return(IT_0306);
}

