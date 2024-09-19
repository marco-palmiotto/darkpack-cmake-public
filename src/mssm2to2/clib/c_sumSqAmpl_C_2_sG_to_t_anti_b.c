#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sG_to_t_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sG_to_t_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = m_b*m_sG*s_13;
    const ccomplex_t IT_0001 = pow(m_t, 2);
    const ccomplex_t IT_0002 = pow(m_C_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = m_b*conj(U_d2)*V_tb*e_em*IT_0010*U_sb_10;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0013 + (-2)*IT_0015);
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0004*IT_0017;
    const ccomplex_t IT_0019 = IT_0003*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_01);
    const ccomplex_t IT_0023 = m_b*conj(U_d2)*V_tb*e_em*IT_0010*U_sb_11;
    const ccomplex_t IT_0024 = IT_0009*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0025 + (-2)*IT_0027);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = IT_0021*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = -IT_0020 + -IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_10);
    const ccomplex_t IT_0035 = IT_0017*IT_0034;
    const ccomplex_t IT_0036 = IT_0003*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_11);
    const ccomplex_t IT_0039 = IT_0029*IT_0038;
    const ccomplex_t IT_0040 = IT_0021*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0037 + IT_0041;
    const ccomplex_t IT_0043 = pow(m_sG, 2);
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0001 + IT_0043 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0006*IT_0008*IT_0010*conj(U_st_01);
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = IT_0044*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0006*IT_0008*IT_0010*conj(U_st_00);
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0001 + IT_0043 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0050 + IT_0057;
    const ccomplex_t IT_0059 = sin(beta);
    const ccomplex_t IT_0060 = cpow(IT_0059, -1);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0008*IT_0010*IT_0060*U_sb_00;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = IT_0034*IT_0062;
    const ccomplex_t IT_0064 = IT_0003*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0008*IT_0010*IT_0060*U_sb_01;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0038*IT_0067;
    const ccomplex_t IT_0069 = IT_0021*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = -IT_0065 + -IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0073 = IT_0008*IT_0060;
    const ccomplex_t IT_0074 = m_t*V_tb*V_u2*e_em*IT_0010*conj(U_st_10);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = V_tb*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0080 = IT_0008*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = IT_0078 + IT_0083;
    const ccomplex_t IT_0085 = IT_0072*IT_0084;
    const ccomplex_t IT_0086 = IT_0055*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0089 = m_t*V_tb*V_u2*e_em*IT_0010*conj(U_st_11);
    const ccomplex_t IT_0090 = IT_0073*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0095 = IT_0008*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = IT_0093 + IT_0098;
    const ccomplex_t IT_0100 = IT_0088*IT_0099;
    const ccomplex_t IT_0101 = IT_0044*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = -IT_0087 + -IT_0102;
    const ccomplex_t IT_0104 = 8*IT_0103;
    const ccomplex_t IT_0105 = IT_0004*IT_0062;
    const ccomplex_t IT_0106 = IT_0003*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = IT_0022*IT_0067;
    const ccomplex_t IT_0109 = IT_0021*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = IT_0107 + IT_0110;
    const ccomplex_t IT_0112 = 16*IT_0111;
    const ccomplex_t IT_0113 = IT_0104 + IT_0112;
    const ccomplex_t IT_0114 = 8*conj(IT_0103);
    const ccomplex_t IT_0115 = 16*conj(IT_0111);
    const ccomplex_t IT_0116 = IT_0114 + IT_0115;
    const ccomplex_t IT_0117 = 8*conj(IT_0033);
    const ccomplex_t IT_0118 = s_13*s_24;
    const ccomplex_t IT_0119 = s_14*s_23;
    const ccomplex_t IT_0120 = s_12*s_34;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = IT_0118 + IT_0119 + IT_0121;
    const ccomplex_t IT_0123 = 8*IT_0111;
    const ccomplex_t IT_0124 = 8*conj(IT_0111);
    const ccomplex_t IT_0125 = m_b*m_t*s_12;
    const ccomplex_t IT_0126 = IT_0047*IT_0088;
    const ccomplex_t IT_0127 = IT_0044*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0052*IT_0072;
    const ccomplex_t IT_0130 = IT_0055*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = -IT_0128 + -IT_0131;
    const ccomplex_t IT_0133 = IT_0053*IT_0084;
    const ccomplex_t IT_0134 = IT_0055*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0045*IT_0099;
    const ccomplex_t IT_0137 = IT_0044*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = IT_0135 + IT_0138;
    const ccomplex_t IT_0140 = 8*IT_0071;
    const ccomplex_t IT_0141 = 8*conj(IT_0071);
    const ccomplex_t IT_0142 = m_sG*s_34*m_C_2;
    const ccomplex_t IT_0143 = 8*IT_0139;
    const ccomplex_t IT_0144 = 8*conj(IT_0139);
    const ccomplex_t IT_0145 = m_t*m_sG*s_14;
    const ccomplex_t IT_0146 = 16*IT_0103;
    const ccomplex_t IT_0147 = IT_0123 + IT_0146;
    const ccomplex_t IT_0148 = 16*conj(IT_0103);
    const ccomplex_t IT_0149 = IT_0124 + IT_0148;
    const ccomplex_t IT_0150 = 16*IT_0132;
    const ccomplex_t IT_0151 = 16*IT_0071;
    const ccomplex_t IT_0152 = m_t*s_24*m_C_2;
    const ccomplex_t IT_0153 = 16*conj(IT_0071);
    const ccomplex_t IT_0154 = 16*IT_0139;
    const ccomplex_t IT_0155 = m_b*m_t*m_sG*m_C_2;
    const ccomplex_t IT_0156 = 16*conj(IT_0139);
    const ccomplex_t IT_0157 = m_b*s_23*m_C_2;
    const ccomplex_t IT_0158 = IT_0000*(16*conj(IT_0033)*IT_0042 + IT_0033*(16
      *conj(IT_0042) + 8*conj(IT_0058)) + conj(IT_0071)*IT_0113 + IT_0071
      *IT_0116 + IT_0058*IT_0117) + 8*IT_0122*(conj(IT_0042)*IT_0058 + IT_0042
      *conj(IT_0058) + 0.125*conj(IT_0103)*IT_0123 + 0.125*IT_0103*IT_0124) + 8
      *IT_0125*(0.125*IT_0117*IT_0132 + IT_0033*conj(IT_0132) + 0.125*conj
      (IT_0139)*IT_0140 + 0.125*IT_0139*IT_0141) + IT_0142*(conj(IT_0132)
      *IT_0140 + IT_0132*IT_0141 + conj(IT_0033)*IT_0143 + IT_0033*IT_0144) + 8
      *IT_0145*(conj(IT_0042)*IT_0132 + (IT_0042 + 2*IT_0058)*conj(IT_0132) +
       0.125*conj(IT_0139)*IT_0147 + 0.125*IT_0139*IT_0149 + 0.125*conj(IT_0058)
      *IT_0150) + IT_0118*(16*IT_0033*conj(IT_0033) + 16*IT_0042*conj(IT_0042) +
       conj(IT_0111)*IT_0112 + conj(IT_0071)*IT_0151) + IT_0152*(conj(IT_0033)
      *IT_0113 + IT_0033*IT_0116 + conj(IT_0058)*IT_0140 + IT_0058*IT_0141 +
       conj(IT_0042)*IT_0151 + IT_0042*IT_0153) + IT_0119*(16*IT_0058*conj
      (IT_0058) + conj(IT_0103)*IT_0146 + conj(IT_0132)*IT_0150 + conj(IT_0139)
      *IT_0154) + IT_0155*(conj(IT_0042)*IT_0113 + IT_0042*IT_0116 + conj
      (IT_0058)*IT_0147 + IT_0058*IT_0149 + conj(IT_0033)*IT_0151 + IT_0033
      *IT_0153 + conj(IT_0132)*IT_0154 + IT_0132*IT_0156) + (conj(IT_0042)
      *IT_0143 + IT_0042*IT_0144 + conj(IT_0132)*IT_0147 + IT_0132*IT_0149 +
       conj(IT_0058)*IT_0154 + IT_0058*IT_0156)*IT_0157;
    return create_ccomplex_return(IT_0158);
}

