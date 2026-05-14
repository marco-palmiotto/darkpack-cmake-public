#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_sG_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_sG_to_anti_b_b(
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
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
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
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = m_b*conj(N_d3)*e_em*IT_0014*U_sb_10;
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0019 = IT_0012*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = cos(theta_W);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0017 + -IT_0020 +
       0.333333333333333*IT_0025);
    const ccomplex_t IT_0027 = 3*IT_0026;
    const ccomplex_t IT_0028 = 0.166666666666667*IT_0027;
    const ccomplex_t IT_0029 = IT_0008*IT_0028;
    const ccomplex_t IT_0030 = IT_0007*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_11);
    const ccomplex_t IT_0034 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0035 = IT_0022*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0038 = IT_0012*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = m_b*conj(N_d3)*e_em*IT_0014*U_sb_11;
    const ccomplex_t IT_0041 = IT_0013*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0036 + (-3)*IT_0039 + 3
      *IT_0042);
    const ccomplex_t IT_0044 = 0.166666666666667*IT_0043;
    const ccomplex_t IT_0045 = IT_0033*IT_0044;
    const ccomplex_t IT_0046 = IT_0032*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = -IT_0031 + -IT_0047;
    const ccomplex_t IT_0049 = pow(m_N_3, 2);
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0005 + IT_0049 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0052 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0053 = IT_0022*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = m_b*conj(N_d3)*e_em*IT_0014*conj(U_sb_00);
    const ccomplex_t IT_0056 = IT_0013*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0054 + 1.5*IT_0057);
    const ccomplex_t IT_0059 = (-0.333333333333333)*IT_0058;
    const ccomplex_t IT_0060 = IT_0051*IT_0059;
    const ccomplex_t IT_0061 = IT_0050*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0005 + IT_0049 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0065 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0066 = IT_0022*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = m_b*conj(N_d3)*e_em*IT_0014*conj(U_sb_01);
    const ccomplex_t IT_0069 = IT_0013*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0067 + 1.5*IT_0070);
    const ccomplex_t IT_0072 = (-0.333333333333333)*IT_0071;
    const ccomplex_t IT_0073 = IT_0064*IT_0072;
    const ccomplex_t IT_0074 = IT_0063*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0062 + IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0078 = m_b*N_d3*e_em*IT_0014*conj(U_sb_10);
    const ccomplex_t IT_0079 = IT_0013*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0082 = IT_0022*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0085 = IT_0012*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0080 + 0.333333333333333
      *IT_0083 + -IT_0086);
    const ccomplex_t IT_0088 = 3*IT_0087;
    const ccomplex_t IT_0089 = 0.166666666666667*IT_0088;
    const ccomplex_t IT_0090 = IT_0077*IT_0089;
    const ccomplex_t IT_0091 = IT_0050*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0094 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0095 = IT_0022*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0098 = IT_0012*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = m_b*N_d3*e_em*IT_0014*conj(U_sb_11);
    const ccomplex_t IT_0101 = IT_0013*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0096 + (-3)*IT_0099 + 3
      *IT_0102);
    const ccomplex_t IT_0104 = 0.166666666666667*IT_0103;
    const ccomplex_t IT_0105 = IT_0093*IT_0104;
    const ccomplex_t IT_0106 = IT_0063*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = IT_0092 + IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0110 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0111 = IT_0022*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = m_b*N_d3*e_em*IT_0014*U_sb_00;
    const ccomplex_t IT_0114 = IT_0013*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0112 + 1.5*IT_0115);
    const ccomplex_t IT_0117 = (-0.333333333333333)*IT_0116;
    const ccomplex_t IT_0118 = IT_0109*IT_0117;
    const ccomplex_t IT_0119 = IT_0007*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_01);
    const ccomplex_t IT_0122 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0123 = IT_0022*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = m_b*N_d3*e_em*IT_0014*U_sb_01;
    const ccomplex_t IT_0126 = IT_0013*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0124 + 1.5*IT_0127);
    const ccomplex_t IT_0129 = (-0.333333333333333)*IT_0128;
    const ccomplex_t IT_0130 = IT_0121*IT_0129;
    const ccomplex_t IT_0131 = IT_0032*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = -IT_0120 + -IT_0132;
    const ccomplex_t IT_0134 = 8*IT_0133;
    const ccomplex_t IT_0135 = 8*conj(IT_0133);
    const ccomplex_t IT_0136 = s_12*IT_0005;
    const ccomplex_t IT_0137 = IT_0051*IT_0089;
    const ccomplex_t IT_0138 = IT_0050*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = IT_0064*IT_0104;
    const ccomplex_t IT_0141 = IT_0063*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = -IT_0139 + -IT_0142;
    const ccomplex_t IT_0144 = IT_0008*IT_0117;
    const ccomplex_t IT_0145 = IT_0007*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0033*IT_0129;
    const ccomplex_t IT_0148 = IT_0032*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0146 + IT_0149;
    const ccomplex_t IT_0151 = 8*IT_0150;
    const ccomplex_t IT_0152 = 8*conj(IT_0150);
    const ccomplex_t IT_0153 = IT_0059*IT_0077;
    const ccomplex_t IT_0154 = IT_0050*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = IT_0072*IT_0093;
    const ccomplex_t IT_0157 = IT_0063*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = -IT_0155 + -IT_0158;
    const ccomplex_t IT_0160 = IT_0028*IT_0109;
    const ccomplex_t IT_0161 = IT_0007*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = IT_0044*IT_0121;
    const ccomplex_t IT_0164 = IT_0032*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = IT_0162 + IT_0165;
    const ccomplex_t IT_0167 = 8*IT_0166;
    const ccomplex_t IT_0168 = 8*conj(IT_0166);
    const ccomplex_t IT_0169 = m_sG*s_34*m_N_3;
    const ccomplex_t IT_0170 = 8*IT_0143;
    const ccomplex_t IT_0171 = 8*conj(IT_0143);
    const ccomplex_t IT_0172 = m_b*s_23*m_N_3;
    const ccomplex_t IT_0173 = 8*IT_0108;
    const ccomplex_t IT_0174 = 16*IT_0133;
    const ccomplex_t IT_0175 = IT_0173 + IT_0174;
    const ccomplex_t IT_0176 = 8*conj(IT_0108);
    const ccomplex_t IT_0177 = 16*conj(IT_0133);
    const ccomplex_t IT_0178 = IT_0176 + IT_0177;
    const ccomplex_t IT_0179 = 16*IT_0150;
    const ccomplex_t IT_0180 = 16*conj(IT_0150);
    const ccomplex_t IT_0181 = 16*IT_0166;
    const ccomplex_t IT_0182 = m_b*m_sG*s_14;
    const ccomplex_t IT_0183 = m_sG*m_N_3*IT_0005;
    const ccomplex_t IT_0184 = 16*IT_0108;
    const ccomplex_t IT_0185 = IT_0134 + IT_0184;
    const ccomplex_t IT_0186 = 16*conj(IT_0108);
    const ccomplex_t IT_0187 = IT_0135 + IT_0186;
    const ccomplex_t IT_0188 = 16*IT_0143;
    const ccomplex_t IT_0189 = 16*conj(IT_0143);
    const ccomplex_t IT_0190 = m_b*s_24*m_N_3;
    const ccomplex_t IT_0191 = 16*IT_0159;
    const ccomplex_t IT_0192 = m_b*m_sG*s_13;
    const ccomplex_t IT_0193 = 8*IT_0004*(conj(IT_0048)*IT_0076 + IT_0048*conj
      (IT_0076) + 0.125*conj(IT_0108)*IT_0134 + 0.125*IT_0108*IT_0135) + IT_0136
      *(conj(IT_0143)*IT_0151 + IT_0143*IT_0152 + conj(IT_0159)*IT_0167 +
       IT_0159*IT_0168) + IT_0169*(IT_0152*IT_0159 + IT_0151*conj(IT_0159) +
       conj(IT_0166)*IT_0170 + IT_0166*IT_0171) + IT_0172*(conj(IT_0076)*IT_0151
       + IT_0076*IT_0152 + conj(IT_0166)*IT_0175 + IT_0166*IT_0178 + conj
      (IT_0048)*IT_0179 + IT_0048*IT_0180) + IT_0000*(16*IT_0048*conj(IT_0048) +
       conj(IT_0133)*IT_0174 + conj(IT_0150)*IT_0179 + conj(IT_0166)*IT_0181) + 
      (16*IT_0048*conj(IT_0166) + conj(IT_0076)*IT_0167 + IT_0076*IT_0168 + conj
      (IT_0150)*IT_0175 + IT_0150*IT_0178 + conj(IT_0048)*IT_0181)*IT_0182 +
       IT_0183*(conj(IT_0048)*IT_0175 + IT_0048*IT_0178 + conj(IT_0166)*IT_0179 
      + IT_0166*IT_0180 + conj(IT_0076)*IT_0185 + IT_0076*IT_0187 + conj(IT_0159
      )*IT_0188 + IT_0159*IT_0189) + (conj(IT_0048)*IT_0170 + IT_0048*IT_0171 +
       conj(IT_0159)*IT_0185 + IT_0159*IT_0187 + conj(IT_0076)*IT_0188 + IT_0076
      *IT_0189)*IT_0190 + IT_0001*(16*IT_0076*conj(IT_0076) + conj(IT_0108)
      *IT_0184 + conj(IT_0143)*IT_0188 + conj(IT_0159)*IT_0191) + 8*(conj
      (IT_0048)*IT_0159 + (IT_0048 + 2*IT_0076)*conj(IT_0159) + 0.125*conj
      (IT_0143)*IT_0185 + 0.125*IT_0143*IT_0187 + 0.125*conj(IT_0076)*IT_0191)
      *IT_0192;
    return create_ccomplex_return(IT_0193);
}

