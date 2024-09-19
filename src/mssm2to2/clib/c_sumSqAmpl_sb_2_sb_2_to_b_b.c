#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_sb_2_to_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_sb_2_to_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_sb_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = cpow((-2)*s_23 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0013 = cos(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = cos(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = sin(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = pow(m_W, -1);
    const ccomplex_t IT_0024 = m_b*N_d1*e_em*IT_0023*U_sb_01;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0017 + 1.5*IT_0026);
    const ccomplex_t IT_0028 = (-0.333333333333333)*IT_0027;
    const ccomplex_t IT_0029 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0030 = IT_0014*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0033 = IT_0021*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = m_b*conj(N_d1)*e_em*IT_0023*U_sb_11;
    const ccomplex_t IT_0036 = IT_0022*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0031 + (-3)*IT_0034 + 3
      *IT_0037);
    const ccomplex_t IT_0039 = 0.166666666666667*IT_0038;
    const ccomplex_t IT_0040 = IT_0028*IT_0039;
    const ccomplex_t IT_0041 = IT_0012*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0044 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0045 = IT_0014*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0048 = IT_0021*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = m_b*conj(N_d2)*e_em*IT_0023*U_sb_11;
    const ccomplex_t IT_0051 = IT_0022*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0046 + (-3)*IT_0049 + 3
      *IT_0052);
    const ccomplex_t IT_0054 = 0.166666666666667*IT_0053;
    const ccomplex_t IT_0055 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0056 = IT_0014*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = m_b*N_d2*e_em*IT_0023*U_sb_01;
    const ccomplex_t IT_0059 = IT_0022*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + 1.5*IT_0060);
    const ccomplex_t IT_0062 = (-0.333333333333333)*IT_0061;
    const ccomplex_t IT_0063 = IT_0054*IT_0062;
    const ccomplex_t IT_0064 = IT_0043*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0067 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0068 = IT_0014*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = m_b*N_d3*e_em*IT_0023*U_sb_01;
    const ccomplex_t IT_0071 = IT_0022*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + 1.5*IT_0072);
    const ccomplex_t IT_0074 = (-0.333333333333333)*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0076 = IT_0014*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0079 = IT_0021*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = m_b*conj(N_d3)*e_em*IT_0023*U_sb_11;
    const ccomplex_t IT_0082 = IT_0022*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0077 + (-3)*IT_0080 + 3
      *IT_0083);
    const ccomplex_t IT_0085 = 0.166666666666667*IT_0084;
    const ccomplex_t IT_0086 = IT_0074*IT_0085;
    const ccomplex_t IT_0087 = IT_0066*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0090 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0091 = IT_0014*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0094 = IT_0021*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = m_b*conj(N_d4)*e_em*IT_0023*U_sb_11;
    const ccomplex_t IT_0097 = IT_0022*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0092 + (-3)*IT_0095 + 3
      *IT_0098);
    const ccomplex_t IT_0100 = 0.166666666666667*IT_0099;
    const ccomplex_t IT_0101 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0102 = IT_0014*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = m_b*N_d4*e_em*IT_0023*U_sb_01;
    const ccomplex_t IT_0105 = IT_0022*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0103 + 1.5*IT_0106);
    const ccomplex_t IT_0108 = (-0.333333333333333)*IT_0107;
    const ccomplex_t IT_0109 = IT_0100*IT_0108;
    const ccomplex_t IT_0110 = IT_0089*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = (-0.5)*IT_0008 + (-0.166666666666667)*IT_0011 +
       IT_0042 + IT_0065 + IT_0088 + IT_0111;
    const ccomplex_t IT_0113 = s_34*IT_0002;
    const ccomplex_t IT_0114 = (-18)*IT_0113;
    const ccomplex_t IT_0115 = s_23*s_24;
    const ccomplex_t IT_0116 = 36*IT_0115;
    const ccomplex_t IT_0117 = IT_0114 + IT_0116;
    const ccomplex_t IT_0118 = IT_0000*IT_0002;
    const ccomplex_t IT_0119 = (-36)*IT_0118;
    const ccomplex_t IT_0120 = 2*IT_0117 + IT_0119;
    const ccomplex_t IT_0121 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0122 = IT_0086*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0125 = IT_0109*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0128 = IT_0040*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0131 = IT_0063*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = 0.166666666666667*IT_0008 + 0.5*IT_0011 + 
      -IT_0123 + -IT_0126 + -IT_0129 + -IT_0132;
    const ccomplex_t IT_0134 = (-6)*IT_0113;
    const ccomplex_t IT_0135 = 2*IT_0134;
    const ccomplex_t IT_0136 = 12*IT_0115;
    const ccomplex_t IT_0137 = 2*IT_0136;
    const ccomplex_t IT_0138 = (-12)*IT_0118;
    const ccomplex_t IT_0139 = IT_0135 + IT_0137 + IT_0138;
    const ccomplex_t IT_0140 = cpow(IT_0005, 2);
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*m_sG*IT_0140;
    const ccomplex_t IT_0142 = IT_0009*IT_0141;
    const ccomplex_t IT_0143 = cpow(IT_0108, 2);
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*m_N_4*IT_0143;
    const ccomplex_t IT_0145 = IT_0124*IT_0144;
    const ccomplex_t IT_0146 = cpow(IT_0028, 2);
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*m_N_1*IT_0146;
    const ccomplex_t IT_0148 = IT_0127*IT_0147;
    const ccomplex_t IT_0149 = cpow(IT_0074, 2);
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*m_N_3*IT_0149;
    const ccomplex_t IT_0151 = IT_0121*IT_0150;
    const ccomplex_t IT_0152 = (-0.5)*IT_0142 + -IT_0145 + -IT_0148 + -IT_0151;
    const ccomplex_t IT_0153 = cpow(IT_0062, 2);
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*m_N_2*IT_0153;
    const ccomplex_t IT_0155 = IT_0043*IT_0154;
    const ccomplex_t IT_0156 = IT_0089*IT_0144;
    const ccomplex_t IT_0157 = IT_0012*IT_0147;
    const ccomplex_t IT_0158 = 0.166666666666667*IT_0142 + -IT_0155 + -IT_0156
       + -IT_0157;
    const ccomplex_t IT_0159 = IT_0003*IT_0141;
    const ccomplex_t IT_0160 = IT_0130*IT_0154;
    const ccomplex_t IT_0161 = m_b*IT_0008;
    const ccomplex_t IT_0162 = m_b*IT_0011;
    const ccomplex_t IT_0163 = m_b*IT_0129;
    const ccomplex_t IT_0164 = m_b*IT_0132;
    const ccomplex_t IT_0165 = m_b*IT_0123;
    const ccomplex_t IT_0166 = m_b*IT_0126;
    const ccomplex_t IT_0167 = (-0.166666666666667)*IT_0159 + IT_0160 +
       0.166666666666667*IT_0161 + (-0.5)*IT_0162 + -IT_0163 + -IT_0164 + 
      -IT_0165 + -IT_0166;
    const ccomplex_t IT_0168 = IT_0066*IT_0150;
    const ccomplex_t IT_0169 = m_b*IT_0042;
    const ccomplex_t IT_0170 = m_b*IT_0065;
    const ccomplex_t IT_0171 = m_b*IT_0088;
    const ccomplex_t IT_0172 = m_b*IT_0111;
    const ccomplex_t IT_0173 = 0.5*IT_0159 + (-0.5)*IT_0161 +
       0.166666666666667*IT_0162 + IT_0168 + -IT_0169 + -IT_0170 + -IT_0171 + 
      -IT_0172;
    const ccomplex_t IT_0174 = cpow(IT_0004, 2);
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*m_sG*IT_0174;
    const ccomplex_t IT_0176 = IT_0009*IT_0175;
    const ccomplex_t IT_0177 = cpow(IT_0100, 2);
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*m_N_4*IT_0177;
    const ccomplex_t IT_0179 = IT_0124*IT_0178;
    const ccomplex_t IT_0180 = cpow(IT_0039, 2);
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*m_N_1*IT_0180;
    const ccomplex_t IT_0182 = IT_0127*IT_0181;
    const ccomplex_t IT_0183 = IT_0003*IT_0175;
    const ccomplex_t IT_0184 = cpow(IT_0085, 2);
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*m_N_3*IT_0184;
    const ccomplex_t IT_0186 = IT_0121*IT_0185;
    const ccomplex_t IT_0187 = (-0.5)*IT_0176 + -IT_0179 + -IT_0182 +
       0.166666666666667*IT_0183 + -IT_0186;
    const ccomplex_t IT_0188 = cpow(IT_0054, 2);
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*m_N_2*IT_0188;
    const ccomplex_t IT_0190 = IT_0043*IT_0189;
    const ccomplex_t IT_0191 = IT_0089*IT_0178;
    const ccomplex_t IT_0192 = 0.166666666666667*IT_0176 + (-0.5)*IT_0183 + 
      -IT_0190 + -IT_0191;
    const ccomplex_t IT_0193 = IT_0130*IT_0189;
    const ccomplex_t IT_0194 = 0.166666666666667*IT_0161 + (-0.5)*IT_0162 + 
      -IT_0163 + -IT_0164 + -IT_0165 + -IT_0166 + IT_0193;
    const ccomplex_t IT_0195 = IT_0012*IT_0181;
    const ccomplex_t IT_0196 = IT_0066*IT_0185;
    const ccomplex_t IT_0197 = (-0.5)*IT_0161 + 0.166666666666667*IT_0162 + 
      -IT_0169 + -IT_0170 + -IT_0171 + -IT_0172 + IT_0195 + IT_0196;
    const ccomplex_t IT_0198 = 6*s_34;
    const ccomplex_t IT_0199 = 18*s_34;
    const ccomplex_t IT_0200 = (-18)*s_34;
    const ccomplex_t IT_0201 = (-6)*s_34;
    const ccomplex_t IT_0202 = m_b*s_24;
    const ccomplex_t IT_0203 = IT_0187 + 3*IT_0192 + -IT_0194 + (-3)*IT_0197;
    const ccomplex_t IT_0204 = 6*IT_0167;
    const ccomplex_t IT_0205 = 18*IT_0173;
    const ccomplex_t IT_0206 = (-18)*IT_0158;
    const ccomplex_t IT_0207 = (-6)*IT_0152;
    const ccomplex_t IT_0208 = (-3)*conj(IT_0197);
    const ccomplex_t IT_0209 = 3*conj(IT_0192);
    const ccomplex_t IT_0210 = -conj(IT_0194);
    const ccomplex_t IT_0211 = conj(IT_0187) + IT_0208 + IT_0209 + IT_0210;
    const ccomplex_t IT_0212 = 6*conj(IT_0167);
    const ccomplex_t IT_0213 = 18*conj(IT_0173);
    const ccomplex_t IT_0214 = (-18)*conj(IT_0158);
    const ccomplex_t IT_0215 = (-6)*conj(IT_0152);
    const ccomplex_t IT_0216 = 6*IT_0173;
    const ccomplex_t IT_0217 = 18*IT_0167;
    const ccomplex_t IT_0218 = (-18)*IT_0152;
    const ccomplex_t IT_0219 = (-6)*IT_0158;
    const ccomplex_t IT_0220 = 6*conj(IT_0173);
    const ccomplex_t IT_0221 = 18*conj(IT_0167);
    const ccomplex_t IT_0222 = (-18)*conj(IT_0152);
    const ccomplex_t IT_0223 = (-6)*conj(IT_0158);
    const ccomplex_t IT_0224 = conj(IT_0112)*(IT_0203 + (-0.166666666666667)
      *IT_0204 + (-0.166666666666667)*IT_0205 + (-0.166666666666667)*IT_0206 + (
      -0.166666666666667)*IT_0207) + IT_0112*(IT_0211 + (-0.166666666666667)
      *IT_0212 + (-0.166666666666667)*IT_0213 + (-0.166666666666667)*IT_0214 + (
      -0.166666666666667)*IT_0215) + (-0.166666666666667)*conj(IT_0133)*((-18)
      *IT_0187 + (-6)*IT_0192 + 18*IT_0194 + 6*IT_0197 + IT_0216 + IT_0217 +
       IT_0218 + IT_0219) + (-0.166666666666667)*IT_0133*((-18)*conj(IT_0187) + 
      (-6)*conj(IT_0192) + 18*conj(IT_0194) + 6*conj(IT_0197) + IT_0220 +
       IT_0221 + IT_0222 + IT_0223);
    const ccomplex_t IT_0225 = 6*IT_0152;
    const ccomplex_t IT_0226 = 18*IT_0158;
    const ccomplex_t IT_0227 = (-18)*IT_0173;
    const ccomplex_t IT_0228 = (-6)*IT_0167;
    const ccomplex_t IT_0229 = 6*conj(IT_0152);
    const ccomplex_t IT_0230 = 18*conj(IT_0158);
    const ccomplex_t IT_0231 = (-18)*conj(IT_0173);
    const ccomplex_t IT_0232 = (-6)*conj(IT_0167);
    const ccomplex_t IT_0233 = 6*IT_0158;
    const ccomplex_t IT_0234 = 18*IT_0152;
    const ccomplex_t IT_0235 = (-18)*IT_0167;
    const ccomplex_t IT_0236 = (-6)*IT_0173;
    const ccomplex_t IT_0237 = 6*conj(IT_0158);
    const ccomplex_t IT_0238 = 18*conj(IT_0152);
    const ccomplex_t IT_0239 = (-18)*conj(IT_0167);
    const ccomplex_t IT_0240 = (-6)*conj(IT_0173);
    const ccomplex_t IT_0241 = conj(IT_0192)*(IT_0204 + IT_0205 + IT_0206 +
       IT_0207) + IT_0192*(IT_0212 + IT_0213 + IT_0214 + IT_0215) + conj(IT_0187
      )*(IT_0216 + IT_0217 + IT_0218 + IT_0219) + IT_0187*(IT_0220 + IT_0221 +
       IT_0222 + IT_0223) + conj(IT_0197)*(IT_0225 + IT_0226 + IT_0227 + IT_0228
      ) + IT_0197*(IT_0229 + IT_0230 + IT_0231 + IT_0232) + conj(IT_0194)*
      (IT_0233 + IT_0234 + IT_0235 + IT_0236) + IT_0194*(IT_0237 + IT_0238 +
       IT_0239 + IT_0240);
    const ccomplex_t IT_0242 = m_b*s_23;
    const ccomplex_t IT_0243 = conj(IT_0112)*(IT_0203 + 0.166666666666667
      *IT_0225 + 0.166666666666667*IT_0226 + 0.166666666666667*IT_0227 +
       0.166666666666667*IT_0228) + IT_0112*(IT_0211 + 0.166666666666667*IT_0229
       + 0.166666666666667*IT_0230 + 0.166666666666667*IT_0231 +
       0.166666666666667*IT_0232) + 0.166666666666667*conj(IT_0133)*(18*IT_0187 
      + 6*IT_0192 + (-18)*IT_0194 + (-6)*IT_0197 + IT_0233 + IT_0234 + IT_0235 +
       IT_0236) + 0.166666666666667*IT_0133*(18*conj(IT_0187) + 6*conj(IT_0192) 
      + (-18)*conj(IT_0194) + (-6)*conj(IT_0197) + IT_0237 + IT_0238 + IT_0239 +
       IT_0240);
    const ccomplex_t IT_0244 = IT_0133*(IT_0120*conj(IT_0133) + conj(IT_0112)
      *IT_0139) + IT_0112*(conj(IT_0112)*IT_0120 + conj(IT_0133)*IT_0139) + 
      (conj(IT_0152)*IT_0158 + IT_0152*conj(IT_0158) + conj(IT_0167)*IT_0173 +
       IT_0167*conj(IT_0173) + conj(IT_0187)*IT_0192 + IT_0187*conj(IT_0192) +
       conj(IT_0194)*IT_0197 + IT_0194*conj(IT_0197))*IT_0198 + (IT_0152*conj
      (IT_0152) + IT_0158*conj(IT_0158) + IT_0167*conj(IT_0167) + IT_0173*conj
      (IT_0173) + IT_0187*conj(IT_0187) + IT_0192*conj(IT_0192) + IT_0194*conj
      (IT_0194) + IT_0197*conj(IT_0197))*IT_0199 + (conj(IT_0152)*IT_0167 +
       IT_0152*conj(IT_0167) + conj(IT_0158)*IT_0173 + IT_0158*conj(IT_0173) +
       conj(IT_0187)*IT_0194 + IT_0187*conj(IT_0194) + conj(IT_0192)*IT_0197 +
       IT_0192*conj(IT_0197))*IT_0200 + (conj(IT_0158)*IT_0167 + IT_0158*conj
      (IT_0167) + conj(IT_0152)*IT_0173 + IT_0152*conj(IT_0173) + conj(IT_0192)
      *IT_0194 + IT_0192*conj(IT_0194) + conj(IT_0187)*IT_0197 + IT_0187*conj
      (IT_0197))*IT_0201 + (-6)*IT_0202*IT_0224 + IT_0000*IT_0241 + 6*IT_0242
      *IT_0243;
    return create_ccomplex_return(IT_0244);
}

