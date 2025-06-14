#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_sG_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_sG_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = s_14*s_23;
    const ccomplex_t IT_0001 = s_13*s_24;
    const ccomplex_t IT_0002 = s_12*s_34;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = IT_0000 + IT_0001 + IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_10);
    const ccomplex_t IT_0006 = cos(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = sin(beta);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = IT_0012*IT_0017;
    const ccomplex_t IT_0019 = pow(m_W, -1);
    const ccomplex_t IT_0020 = m_t*conj(N_u2)*e_em*IT_0019*U_st_10;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0010 + 3*IT_0015 + 3
      *IT_0022);
    const ccomplex_t IT_0024 = 0.166666666666667*IT_0023;
    const ccomplex_t IT_0025 = IT_0005*IT_0024;
    const ccomplex_t IT_0026 = pow(m_t, 2);
    const ccomplex_t IT_0027 = pow(m_sG, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_11);
    const ccomplex_t IT_0032 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0033 = IT_0007*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0036 = IT_0012*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = m_t*conj(N_u2)*e_em*IT_0019*U_st_11;
    const ccomplex_t IT_0039 = IT_0018*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0034 + 3*IT_0037 + 3
      *IT_0040);
    const ccomplex_t IT_0042 = 0.166666666666667*IT_0041;
    const ccomplex_t IT_0043 = IT_0031*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = -IT_0030 + -IT_0046;
    const ccomplex_t IT_0048 = conj(N_B2)*e_em*conj(U_st_11);
    const ccomplex_t IT_0049 = IT_0007*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = m_t*conj(N_u2)*e_em*IT_0019*conj(U_st_01);
    const ccomplex_t IT_0052 = IT_0018*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + (-0.75)*IT_0053);
    const ccomplex_t IT_0055 = 0.666666666666667*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = pow(m_N_2, 2);
    const ccomplex_t IT_0059 = cpow((-2)*s_13 + IT_0026 + IT_0058 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_13 + IT_0026 + IT_0058 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0064 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0065 = IT_0007*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = m_t*conj(N_u2)*e_em*IT_0019*conj(U_st_00);
    const ccomplex_t IT_0068 = IT_0018*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + (-0.75)*IT_0069);
    const ccomplex_t IT_0071 = 0.666666666666667*IT_0070;
    const ccomplex_t IT_0072 = IT_0063*IT_0071;
    const ccomplex_t IT_0073 = IT_0062*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0061 + IT_0074;
    const ccomplex_t IT_0076 = N_B2*e_em*conj(U_st_01);
    const ccomplex_t IT_0077 = IT_0007*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = N_W2*e_em*conj(U_st_01);
    const ccomplex_t IT_0080 = IT_0012*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = m_t*N_u2*e_em*IT_0019*conj(U_st_11);
    const ccomplex_t IT_0083 = IT_0018*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0078 + 3*IT_0081 + 3
      *IT_0084);
    const ccomplex_t IT_0086 = 0.166666666666667*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = IT_0059*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0092 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0093 = IT_0007*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0096 = IT_0012*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = m_t*N_u2*e_em*IT_0019*conj(U_st_10);
    const ccomplex_t IT_0099 = IT_0018*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0094 + 3*IT_0097 + 3
      *IT_0100);
    const ccomplex_t IT_0102 = 0.166666666666667*IT_0101;
    const ccomplex_t IT_0103 = IT_0091*IT_0102;
    const ccomplex_t IT_0104 = IT_0062*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0090 + IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_00);
    const ccomplex_t IT_0108 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0109 = IT_0007*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = m_t*N_u2*e_em*IT_0019*U_st_00;
    const ccomplex_t IT_0112 = IT_0018*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + (-0.75)*IT_0113);
    const ccomplex_t IT_0115 = 0.666666666666667*IT_0114;
    const ccomplex_t IT_0116 = IT_0107*IT_0115;
    const ccomplex_t IT_0117 = IT_0028*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_01);
    const ccomplex_t IT_0120 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0121 = IT_0007*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = m_t*N_u2*e_em*IT_0019*U_st_01;
    const ccomplex_t IT_0124 = IT_0018*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0122 + (-0.75)*IT_0125);
    const ccomplex_t IT_0127 = 0.666666666666667*IT_0126;
    const ccomplex_t IT_0128 = IT_0119*IT_0127;
    const ccomplex_t IT_0129 = IT_0044*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = -IT_0118 + -IT_0130;
    const ccomplex_t IT_0132 = 8*IT_0131;
    const ccomplex_t IT_0133 = 8*conj(IT_0131);
    const ccomplex_t IT_0134 = m_sG*s_34*m_N_2;
    const ccomplex_t IT_0135 = IT_0055*IT_0087;
    const ccomplex_t IT_0136 = IT_0059*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = IT_0071*IT_0091;
    const ccomplex_t IT_0139 = IT_0062*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = -IT_0137 + -IT_0140;
    const ccomplex_t IT_0142 = IT_0031*IT_0127;
    const ccomplex_t IT_0143 = IT_0044*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = IT_0005*IT_0115;
    const ccomplex_t IT_0146 = IT_0028*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = IT_0144 + IT_0147;
    const ccomplex_t IT_0149 = 8*IT_0148;
    const ccomplex_t IT_0150 = 8*conj(IT_0148);
    const ccomplex_t IT_0151 = IT_0024*IT_0107;
    const ccomplex_t IT_0152 = IT_0028*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = IT_0042*IT_0119;
    const ccomplex_t IT_0155 = IT_0044*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = IT_0153 + IT_0156;
    const ccomplex_t IT_0158 = IT_0056*IT_0086;
    const ccomplex_t IT_0159 = IT_0059*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = IT_0063*IT_0102;
    const ccomplex_t IT_0162 = IT_0062*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = -IT_0160 + -IT_0163;
    const ccomplex_t IT_0165 = 8*IT_0164;
    const ccomplex_t IT_0166 = 8*conj(IT_0164);
    const ccomplex_t IT_0167 = s_12*IT_0026;
    const ccomplex_t IT_0168 = 8*IT_0157;
    const ccomplex_t IT_0169 = 8*conj(IT_0157);
    const ccomplex_t IT_0170 = m_t*s_23*m_N_2;
    const ccomplex_t IT_0171 = 8*IT_0106;
    const ccomplex_t IT_0172 = 16*IT_0131;
    const ccomplex_t IT_0173 = IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = 8*conj(IT_0106);
    const ccomplex_t IT_0175 = 16*conj(IT_0131);
    const ccomplex_t IT_0176 = IT_0174 + IT_0175;
    const ccomplex_t IT_0177 = 16*IT_0148;
    const ccomplex_t IT_0178 = 16*conj(IT_0148);
    const ccomplex_t IT_0179 = m_sG*m_N_2*IT_0026;
    const ccomplex_t IT_0180 = 16*IT_0106;
    const ccomplex_t IT_0181 = IT_0132 + IT_0180;
    const ccomplex_t IT_0182 = 16*conj(IT_0106);
    const ccomplex_t IT_0183 = IT_0133 + IT_0182;
    const ccomplex_t IT_0184 = 16*IT_0164;
    const ccomplex_t IT_0185 = 16*conj(IT_0164);
    const ccomplex_t IT_0186 = m_t*s_24*m_N_2;
    const ccomplex_t IT_0187 = 16*IT_0157;
    const ccomplex_t IT_0188 = m_t*m_sG*s_14;
    const ccomplex_t IT_0189 = m_t*m_sG*s_13;
    const ccomplex_t IT_0190 = 16*IT_0141;
    const ccomplex_t IT_0191 = 8*IT_0004*(conj(IT_0047)*IT_0075 + IT_0047*conj
      (IT_0075) + 0.125*conj(IT_0106)*IT_0132 + 0.125*IT_0106*IT_0133) + IT_0134
      *(conj(IT_0141)*IT_0149 + IT_0141*IT_0150 + conj(IT_0157)*IT_0165 +
       IT_0157*IT_0166) + IT_0167*(IT_0150*IT_0164 + IT_0149*conj(IT_0164) +
       conj(IT_0141)*IT_0168 + IT_0141*IT_0169) + IT_0170*(conj(IT_0075)*IT_0149
       + IT_0075*IT_0150 + conj(IT_0157)*IT_0173 + IT_0157*IT_0176 + conj
      (IT_0047)*IT_0177 + IT_0047*IT_0178) + IT_0179*(conj(IT_0047)*IT_0173 +
       IT_0047*IT_0176 + conj(IT_0157)*IT_0177 + IT_0157*IT_0178 + conj(IT_0075)
      *IT_0181 + IT_0075*IT_0183 + conj(IT_0141)*IT_0184 + IT_0141*IT_0185) + 
      (conj(IT_0047)*IT_0165 + IT_0047*IT_0166 + conj(IT_0141)*IT_0181 + IT_0141
      *IT_0183 + conj(IT_0075)*IT_0184 + IT_0075*IT_0185)*IT_0186 + IT_0000*(16
      *IT_0047*conj(IT_0047) + conj(IT_0131)*IT_0172 + conj(IT_0148)*IT_0177 +
       conj(IT_0157)*IT_0187) + (16*IT_0047*conj(IT_0157) + conj(IT_0075)
      *IT_0168 + IT_0075*IT_0169 + conj(IT_0148)*IT_0173 + IT_0148*IT_0176 +
       conj(IT_0047)*IT_0187)*IT_0188 + 8*IT_0189*(conj(IT_0047)*IT_0141 + 
      (IT_0047 + 2*IT_0075)*conj(IT_0141) + 0.125*conj(IT_0164)*IT_0181 + 0.125
      *IT_0164*IT_0183 + 0.125*conj(IT_0075)*IT_0190) + IT_0001*(16*IT_0075*conj
      (IT_0075) + conj(IT_0106)*IT_0180 + conj(IT_0164)*IT_0184 + conj(IT_0141)
      *IT_0190);
    return create_ccomplex_return(IT_0191);
}

