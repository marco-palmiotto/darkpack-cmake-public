#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_sG_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_sG_to_anti_b_b(
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
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = s_14*s_23;
    const ccomplex_t IT_0001 = s_13*s_24;
    const ccomplex_t IT_0002 = s_12*s_34;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = IT_0000 + IT_0001 + IT_0003;
    const ccomplex_t IT_0005 = pow(m_b, 2);
    const ccomplex_t IT_0006 = pow(m_sG, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_10);
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = cos(beta);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = IT_0015*IT_0020;
    const ccomplex_t IT_0022 = pow(m_W, -1);
    const ccomplex_t IT_0023 = m_b*conj(N_d1)*e_em*IT_0022*U_sb_10;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0013 + (-3)*IT_0018 + 3
      *IT_0025);
    const ccomplex_t IT_0027 = (-0.333333333333333)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = 3*IT_0028;
    const ccomplex_t IT_0030 = 0.166666666666667*IT_0029;
    const ccomplex_t IT_0031 = IT_0008*IT_0030;
    const ccomplex_t IT_0032 = IT_0007*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_11);
    const ccomplex_t IT_0036 = m_b*conj(N_d1)*e_em*IT_0022*U_sb_11;
    const ccomplex_t IT_0037 = IT_0021*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0040 = IT_0015*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0043 = IT_0010*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0038 + -IT_0041 +
       0.333333333333333*IT_0044);
    const ccomplex_t IT_0046 = 3*IT_0045;
    const ccomplex_t IT_0047 = 0.166666666666667*IT_0046;
    const ccomplex_t IT_0048 = IT_0035*IT_0047;
    const ccomplex_t IT_0049 = IT_0034*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = -IT_0033 + -IT_0050;
    const ccomplex_t IT_0052 = pow(m_N_1, 2);
    const ccomplex_t IT_0053 = cpow((-2)*s_13 + IT_0005 + IT_0052 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0055 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0056 = IT_0010*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = m_b*conj(N_d1)*e_em*IT_0022*conj(U_sb_00);
    const ccomplex_t IT_0059 = IT_0021*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + 1.5*IT_0060);
    const ccomplex_t IT_0062 = (-0.333333333333333)*IT_0061;
    const ccomplex_t IT_0063 = IT_0054*IT_0062;
    const ccomplex_t IT_0064 = IT_0053*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_13 + IT_0005 + IT_0052 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0068 = conj(N_B1)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0069 = IT_0010*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = m_b*conj(N_d1)*e_em*IT_0022*conj(U_sb_01);
    const ccomplex_t IT_0072 = IT_0021*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0070 + 1.5*IT_0073);
    const ccomplex_t IT_0075 = (-0.333333333333333)*IT_0074;
    const ccomplex_t IT_0076 = IT_0067*IT_0075;
    const ccomplex_t IT_0077 = IT_0066*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0065 + IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0081 = m_b*N_d1*e_em*IT_0022*conj(U_sb_10);
    const ccomplex_t IT_0082 = IT_0021*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0085 = IT_0010*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0088 = IT_0015*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0083 + 0.333333333333333
      *IT_0086 + -IT_0089);
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = (-3)*IT_0091;
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = IT_0080*IT_0093;
    const ccomplex_t IT_0095 = IT_0053*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0098 = N_W1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0099 = IT_0015*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = m_b*N_d1*e_em*IT_0022*conj(U_sb_11);
    const ccomplex_t IT_0102 = IT_0021*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = N_B1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0105 = IT_0010*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0100 + -IT_0103 + (
      -0.333333333333333)*IT_0106);
    const ccomplex_t IT_0108 = (-3)*IT_0107;
    const ccomplex_t IT_0109 = 0.166666666666667*IT_0108;
    const ccomplex_t IT_0110 = IT_0097*IT_0109;
    const ccomplex_t IT_0111 = IT_0066*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = IT_0096 + IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0115 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0116 = IT_0010*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = m_b*N_d1*e_em*IT_0022*U_sb_00;
    const ccomplex_t IT_0119 = IT_0021*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*(IT_0117 + 1.5*IT_0120);
    const ccomplex_t IT_0122 = (-0.333333333333333)*IT_0121;
    const ccomplex_t IT_0123 = IT_0114*IT_0122;
    const ccomplex_t IT_0124 = IT_0007*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_01);
    const ccomplex_t IT_0127 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0128 = IT_0010*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = m_b*N_d1*e_em*IT_0022*U_sb_01;
    const ccomplex_t IT_0131 = IT_0021*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0129 + 1.5*IT_0132);
    const ccomplex_t IT_0134 = (-0.333333333333333)*IT_0133;
    const ccomplex_t IT_0135 = IT_0126*IT_0134;
    const ccomplex_t IT_0136 = IT_0034*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = -IT_0125 + -IT_0137;
    const ccomplex_t IT_0139 = 8*IT_0138;
    const ccomplex_t IT_0140 = 8*conj(IT_0138);
    const ccomplex_t IT_0141 = s_12*IT_0005;
    const ccomplex_t IT_0142 = IT_0054*IT_0093;
    const ccomplex_t IT_0143 = IT_0053*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = IT_0067*IT_0109;
    const ccomplex_t IT_0146 = IT_0066*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = -IT_0144 + -IT_0147;
    const ccomplex_t IT_0149 = IT_0008*IT_0122;
    const ccomplex_t IT_0150 = IT_0007*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = IT_0035*IT_0134;
    const ccomplex_t IT_0153 = IT_0034*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0151 + IT_0154;
    const ccomplex_t IT_0156 = 8*IT_0155;
    const ccomplex_t IT_0157 = 8*conj(IT_0155);
    const ccomplex_t IT_0158 = IT_0075*IT_0097;
    const ccomplex_t IT_0159 = IT_0066*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = IT_0062*IT_0080;
    const ccomplex_t IT_0162 = IT_0053*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = -IT_0160 + -IT_0163;
    const ccomplex_t IT_0165 = IT_0030*IT_0114;
    const ccomplex_t IT_0166 = IT_0007*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = IT_0047*IT_0126;
    const ccomplex_t IT_0169 = IT_0034*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = IT_0167 + IT_0170;
    const ccomplex_t IT_0172 = 8*IT_0171;
    const ccomplex_t IT_0173 = 8*conj(IT_0171);
    const ccomplex_t IT_0174 = m_sG*s_34*m_N_1;
    const ccomplex_t IT_0175 = 8*IT_0148;
    const ccomplex_t IT_0176 = 8*conj(IT_0148);
    const ccomplex_t IT_0177 = m_b*s_23*m_N_1;
    const ccomplex_t IT_0178 = 8*IT_0113;
    const ccomplex_t IT_0179 = 16*IT_0138;
    const ccomplex_t IT_0180 = IT_0178 + IT_0179;
    const ccomplex_t IT_0181 = 8*conj(IT_0113);
    const ccomplex_t IT_0182 = 16*conj(IT_0138);
    const ccomplex_t IT_0183 = IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = 16*IT_0155;
    const ccomplex_t IT_0185 = 16*conj(IT_0155);
    const ccomplex_t IT_0186 = 16*IT_0171;
    const ccomplex_t IT_0187 = m_b*m_sG*s_14;
    const ccomplex_t IT_0188 = m_sG*m_N_1*IT_0005;
    const ccomplex_t IT_0189 = 16*IT_0113;
    const ccomplex_t IT_0190 = IT_0139 + IT_0189;
    const ccomplex_t IT_0191 = 16*conj(IT_0113);
    const ccomplex_t IT_0192 = IT_0140 + IT_0191;
    const ccomplex_t IT_0193 = 16*IT_0148;
    const ccomplex_t IT_0194 = 16*conj(IT_0148);
    const ccomplex_t IT_0195 = m_b*s_24*m_N_1;
    const ccomplex_t IT_0196 = 16*IT_0164;
    const ccomplex_t IT_0197 = m_b*m_sG*s_13;
    const ccomplex_t IT_0198 = 8*IT_0004*(conj(IT_0051)*IT_0079 + IT_0051*conj
      (IT_0079) + 0.125*conj(IT_0113)*IT_0139 + 0.125*IT_0113*IT_0140) + IT_0141
      *(conj(IT_0148)*IT_0156 + IT_0148*IT_0157 + conj(IT_0164)*IT_0172 +
       IT_0164*IT_0173) + IT_0174*(IT_0157*IT_0164 + IT_0156*conj(IT_0164) +
       conj(IT_0171)*IT_0175 + IT_0171*IT_0176) + IT_0177*(conj(IT_0079)*IT_0156
       + IT_0079*IT_0157 + conj(IT_0171)*IT_0180 + IT_0171*IT_0183 + conj
      (IT_0051)*IT_0184 + IT_0051*IT_0185) + IT_0000*(16*IT_0051*conj(IT_0051) +
       conj(IT_0138)*IT_0179 + conj(IT_0155)*IT_0184 + conj(IT_0171)*IT_0186) + 
      (16*IT_0051*conj(IT_0171) + conj(IT_0079)*IT_0172 + IT_0079*IT_0173 + conj
      (IT_0155)*IT_0180 + IT_0155*IT_0183 + conj(IT_0051)*IT_0186)*IT_0187 +
       IT_0188*(conj(IT_0051)*IT_0180 + IT_0051*IT_0183 + conj(IT_0171)*IT_0184 
      + IT_0171*IT_0185 + conj(IT_0079)*IT_0190 + IT_0079*IT_0192 + conj(IT_0164
      )*IT_0193 + IT_0164*IT_0194) + (conj(IT_0051)*IT_0175 + IT_0051*IT_0176 +
       conj(IT_0164)*IT_0190 + IT_0164*IT_0192 + conj(IT_0079)*IT_0193 + IT_0079
      *IT_0194)*IT_0195 + IT_0001*(16*IT_0079*conj(IT_0079) + conj(IT_0113)
      *IT_0189 + conj(IT_0148)*IT_0193 + conj(IT_0164)*IT_0196) + 8*(conj
      (IT_0051)*IT_0164 + (IT_0051 + 2*IT_0079)*conj(IT_0164) + 0.125*conj
      (IT_0148)*IT_0190 + 0.125*IT_0148*IT_0192 + 0.125*conj(IT_0079)*IT_0196)
      *IT_0197;
    return create_ccomplex_return(IT_0198);
}

