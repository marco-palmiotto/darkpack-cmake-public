#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_sc_R_to_c_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_sc_R_to_c_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0008;
    const ccomplex_t IT_0010 = 0.666666666666667*IT_0009;
    const ccomplex_t IT_0011 = IT_0006*IT_0010;
    const ccomplex_t IT_0012 = pow(m_c, 2);
    const ccomplex_t IT_0013 = pow(m_sc_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0008;
    const ccomplex_t IT_0020 = 0.666666666666667*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0008;
    const ccomplex_t IT_0028 = 0.666666666666667*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0008;
    const ccomplex_t IT_0036 = 0.666666666666667*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = -IT_0016 + -IT_0024 + -IT_0032 + -IT_0040;
    const ccomplex_t IT_0042 = sin(beta);
    const ccomplex_t IT_0043 = cpow(IT_0042, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0000*IT_0004*IT_0043;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = (-0.333333333333333)*IT_0027;
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = IT_0030*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0000*IT_0004*IT_0043;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = (-0.333333333333333)*IT_0019;
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = IT_0022*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0000*IT_0004*IT_0043;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = (-0.333333333333333)*IT_0035;
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = IT_0038*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0000*IT_0004*IT_0043;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = (-0.333333333333333)*IT_0009;
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = IT_0014*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0049 + -IT_0055 + -IT_0061 + -IT_0067;
    const ccomplex_t IT_0069 = pow(m_sd_R, 2);
    const ccomplex_t IT_0070 = m_c*m_d*IT_0069;
    const ccomplex_t IT_0071 = s_34*IT_0069;
    const ccomplex_t IT_0072 = (-18)*IT_0071;
    const ccomplex_t IT_0073 = s_13*s_14;
    const ccomplex_t IT_0074 = 36*IT_0073;
    const ccomplex_t IT_0075 = IT_0072 + IT_0074;
    const ccomplex_t IT_0076 = m_d*s_13;
    const ccomplex_t IT_0077 = m_d*IT_0049;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0080 = IT_0034*IT_0057;
    const ccomplex_t IT_0081 = IT_0038*IT_0079*IT_0080;
    const ccomplex_t IT_0082 = m_d*IT_0055;
    const ccomplex_t IT_0083 = m_d*IT_0061;
    const ccomplex_t IT_0084 = IT_0081 + -IT_0082 + -IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0086 = cpow((-2)*s_13 + IT_0012 + IT_0069 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_cd*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *conj(V_u2)*e_em*IT_0000*IT_0004*IT_0043;
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = IT_0085*IT_0086*IT_0091;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = m_d*IT_0067;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0097 = IT_0006*IT_0063;
    const ccomplex_t IT_0098 = IT_0014*IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0100 = IT_0018*IT_0051;
    const ccomplex_t IT_0101 = IT_0022*IT_0099*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0103 = IT_0026*IT_0045;
    const ccomplex_t IT_0104 = IT_0030*IT_0102*IT_0103;
    const ccomplex_t IT_0105 = -IT_0098 + -IT_0101 + -IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0107 = cpow((-2)*s_13 + IT_0012 + IT_0069 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_cd*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *conj(V_u1)*e_em*IT_0000*IT_0004*IT_0043;
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0109*IT_0111;
    const ccomplex_t IT_0113 = IT_0106*IT_0107*IT_0112;
    const ccomplex_t IT_0114 = pow(m_sG, 2);
    const ccomplex_t IT_0115 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -IT_0114 +
       reg_prop, -1);
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0117 = cpow(IT_0116, 2);
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_sG*IT_0117;
    const ccomplex_t IT_0119 = IT_0115*IT_0118;
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = 6*IT_0120;
    const ccomplex_t IT_0122 = IT_0010*IT_0064;
    const ccomplex_t IT_0123 = IT_0014*IT_0096*IT_0122;
    const ccomplex_t IT_0124 = 0.166666666666667*IT_0119 + -IT_0123;
    const ccomplex_t IT_0125 = 18*IT_0124;
    const ccomplex_t IT_0126 = IT_0036*IT_0058;
    const ccomplex_t IT_0127 = IT_0038*IT_0079*IT_0126;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = 18*IT_0128;
    const ccomplex_t IT_0130 = IT_0028*IT_0046;
    const ccomplex_t IT_0131 = IT_0030*IT_0102*IT_0130;
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = 18*IT_0132;
    const ccomplex_t IT_0134 = m_d*IT_0016;
    const ccomplex_t IT_0135 = -IT_0134;
    const ccomplex_t IT_0136 = (-18)*IT_0135;
    const ccomplex_t IT_0137 = IT_0020*IT_0052;
    const ccomplex_t IT_0138 = IT_0022*IT_0099*IT_0137;
    const ccomplex_t IT_0139 = m_d*IT_0032;
    const ccomplex_t IT_0140 = m_d*IT_0024;
    const ccomplex_t IT_0141 = IT_0138 + -IT_0139 + -IT_0140;
    const ccomplex_t IT_0142 = (-18)*IT_0141;
    const ccomplex_t IT_0143 = m_d*IT_0040;
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = (-18)*IT_0144;
    const ccomplex_t IT_0146 = IT_0121 + IT_0125 + IT_0129 + IT_0133 + IT_0136
       + IT_0142 + IT_0145;
    const ccomplex_t IT_0147 = 18*conj(IT_0124);
    const ccomplex_t IT_0148 = 18*conj(IT_0128);
    const ccomplex_t IT_0149 = 18*conj(IT_0132);
    const ccomplex_t IT_0150 = (-18)*conj(IT_0135);
    const ccomplex_t IT_0151 = (-18)*conj(IT_0141);
    const ccomplex_t IT_0152 = (-18)*conj(IT_0144);
    const ccomplex_t IT_0153 = (-6)*IT_0120;
    const ccomplex_t IT_0154 = IT_0147 + IT_0148 + IT_0149 + IT_0150 + IT_0151
       + IT_0152 + IT_0153;
    const ccomplex_t IT_0155 = conj(IT_0068)*((-18)*IT_0078 + (-18)*IT_0084 +
       6*IT_0093 + (-18)*IT_0095 + 18*IT_0105 + (-6)*IT_0113) + IT_0068*((-18)
      *conj(IT_0078) + (-18)*conj(IT_0084) + 6*conj(IT_0093) + (-18)*conj
      (IT_0095) + 18*conj(IT_0105) + (-6)*conj(IT_0113)) + conj(IT_0041)*IT_0146
       + IT_0041*IT_0154;
    const ccomplex_t IT_0156 = IT_0078 + IT_0084 + IT_0095;
    const ccomplex_t IT_0157 = conj(IT_0078) + conj(IT_0084) + conj(IT_0095);
    const ccomplex_t IT_0158 = 6*s_34;
    const ccomplex_t IT_0159 = conj(IT_0141) + conj(IT_0144);
    const ccomplex_t IT_0160 = conj(IT_0128) + conj(IT_0132);
    const ccomplex_t IT_0161 = conj(IT_0084) + conj(IT_0095);
    const ccomplex_t IT_0162 = 18*s_34;
    const ccomplex_t IT_0163 = IT_0078*IT_0162;
    const ccomplex_t IT_0164 = conj(IT_0078)*IT_0162;
    const ccomplex_t IT_0165 = IT_0084*IT_0162;
    const ccomplex_t IT_0166 = cpow(IT_0120, 2);
    const ccomplex_t IT_0167 = IT_0128 + IT_0132;
    const ccomplex_t IT_0168 = (-18)*s_34;
    const ccomplex_t IT_0169 = IT_0124*IT_0162;
    const ccomplex_t IT_0170 = IT_0135*IT_0168;
    const ccomplex_t IT_0171 = conj(IT_0135)*IT_0168;
    const ccomplex_t IT_0172 = IT_0141*IT_0162;
    const ccomplex_t IT_0173 = IT_0124*IT_0168;
    const ccomplex_t IT_0174 = conj(IT_0135)*IT_0162;
    const ccomplex_t IT_0175 = conj(IT_0124)*IT_0168;
    const ccomplex_t IT_0176 = IT_0128*IT_0162;
    const ccomplex_t IT_0177 = IT_0141*IT_0168;
    const ccomplex_t IT_0178 = conj(IT_0124)*IT_0162;
    const ccomplex_t IT_0179 = conj(IT_0141)*IT_0168;
    const ccomplex_t IT_0180 = m_c*s_14;
    const ccomplex_t IT_0181 = 18*IT_0135;
    const ccomplex_t IT_0182 = 18*IT_0141;
    const ccomplex_t IT_0183 = 18*IT_0144;
    const ccomplex_t IT_0184 = (-18)*IT_0124;
    const ccomplex_t IT_0185 = (-18)*IT_0128;
    const ccomplex_t IT_0186 = (-18)*IT_0132;
    const ccomplex_t IT_0187 = IT_0153 + IT_0181 + IT_0182 + IT_0183 + IT_0184
       + IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = 18*conj(IT_0135);
    const ccomplex_t IT_0189 = 18*conj(IT_0141);
    const ccomplex_t IT_0190 = 18*conj(IT_0144);
    const ccomplex_t IT_0191 = (-18)*conj(IT_0124);
    const ccomplex_t IT_0192 = (-18)*conj(IT_0128);
    const ccomplex_t IT_0193 = (-18)*conj(IT_0132);
    const ccomplex_t IT_0194 = conj(IT_0041)*(18*IT_0078 + 18*IT_0084 + (-6)
      *IT_0093 + 18*IT_0095 + (-18)*IT_0105 + 6*IT_0113) + IT_0041*(18*conj
      (IT_0078) + 18*conj(IT_0084) + (-6)*conj(IT_0093) + 18*conj(IT_0095) + (
      -18)*conj(IT_0105) + 6*conj(IT_0113)) + conj(IT_0068)*IT_0187 + IT_0068*
      (IT_0121 + IT_0188 + IT_0189 + IT_0190 + IT_0191 + IT_0192 + IT_0193);
    const ccomplex_t IT_0195 = m_c*m_d;
    const ccomplex_t IT_0196 = -conj(IT_0132);
    const ccomplex_t IT_0197 = -conj(IT_0128);
    const ccomplex_t IT_0198 = -conj(IT_0124);
    const ccomplex_t IT_0199 = conj(IT_0135) + conj(IT_0141) + conj(IT_0144) +
       IT_0196 + IT_0197 + IT_0198;
    const ccomplex_t IT_0200 = 18*IT_0120;
    const ccomplex_t IT_0201 = -IT_0132;
    const ccomplex_t IT_0202 = -IT_0128;
    const ccomplex_t IT_0203 = -IT_0124 + IT_0135 + IT_0141 + IT_0144 +
       IT_0201 + IT_0202;
    const ccomplex_t IT_0204 = (-18)*IT_0120;
    const ccomplex_t IT_0205 = IT_0154*IT_0156 + IT_0146*IT_0157 + conj
      (IT_0105)*IT_0187 + IT_0105*(IT_0121 + IT_0188 + IT_0189 + IT_0190 +
       IT_0191 + IT_0192 + IT_0193) + 6*IT_0093*(IT_0199 + 0.166666666666667
      *IT_0200) + (-6)*conj(IT_0113)*((-0.166666666666667)*IT_0200 + IT_0203) +
       6*conj(IT_0093)*(IT_0203 + 0.166666666666667*IT_0204) + (-6)*IT_0113*
      (IT_0199 + (-0.166666666666667)*IT_0204);
    const ccomplex_t IT_0206 = (-6)*s_34;
    const ccomplex_t IT_0207 = conj(IT_0041)*((-18)*IT_0068*IT_0070 + IT_0041
      *IT_0075) + conj(IT_0068)*((-18)*IT_0041*IT_0070 + IT_0068*IT_0075) +
       IT_0076*IT_0155 + (conj(IT_0093)*IT_0105 + IT_0093*conj(IT_0105) +
       IT_0120*(conj(IT_0124) + conj(IT_0128) + conj(IT_0132) + IT_0135 +
       IT_0141 + IT_0144) + conj(IT_0113)*IT_0156 + IT_0113*IT_0157)*IT_0158 + 
      (IT_0093*conj(IT_0093) + IT_0105*conj(IT_0105) + IT_0113*conj(IT_0113) +
       IT_0135*(conj(IT_0135) + conj(IT_0141) + conj(IT_0144)) + IT_0144*IT_0159
       + IT_0132*IT_0160 + IT_0095*IT_0161)*IT_0162 + IT_0157*IT_0163 + (IT_0084
       + IT_0095)*IT_0164 + IT_0161*IT_0165 + (conj(IT_0093)*IT_0113 + IT_0093
      *conj(IT_0113) + conj(IT_0105)*IT_0156 + IT_0105*IT_0157 + IT_0144*IT_0160
       + IT_0166 + conj(IT_0144)*IT_0167)*IT_0168 + conj(IT_0124)*(IT_0169 +
       IT_0170) + IT_0124*IT_0171 + IT_0159*(IT_0172 + IT_0173) + (IT_0141 +
       IT_0144)*(IT_0174 + IT_0175) + IT_0160*(IT_0169 + IT_0170 + IT_0176 +
       IT_0177) + IT_0167*(IT_0171 + IT_0178 + IT_0179) + IT_0180*IT_0194 +
       IT_0195*IT_0205 + (conj(IT_0105)*IT_0113 + IT_0105*conj(IT_0113) +
       IT_0120*(IT_0124 + IT_0128 + IT_0132 + conj(IT_0135) + conj(IT_0141) +
       conj(IT_0144)) + conj(IT_0093)*IT_0156 + IT_0093*IT_0157)*IT_0206;
    return create_ccomplex_return(IT_0207);
}

