#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_anti_sb_2_to_s_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_anti_sb_2_to_s_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_ss_R = param->m_ss_R;
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
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_ss_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_11);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0011 = pow(m_W, -1);
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0011*IT_0013*IT_0015;
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = N_W1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0015;
    const ccomplex_t IT_0022 = m_b*N_d1*e_em*IT_0011*conj(U_sb_11);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = cos(theta_W);
    const ccomplex_t IT_0026 = cpow(IT_0025, -1);
    const ccomplex_t IT_0027 = N_B1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0020 + -IT_0024 + (
      -0.333333333333333)*IT_0029);
    const ccomplex_t IT_0031 = (-3)*IT_0030;
    const ccomplex_t IT_0032 = 0.166666666666667*IT_0031;
    const ccomplex_t IT_0033 = IT_0017*IT_0032;
    const ccomplex_t IT_0034 = IT_0010*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0011*IT_0013*IT_0015;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = N_B2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0040 = IT_0026*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = N_W2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0043 = IT_0015*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = m_b*N_d2*e_em*IT_0011*conj(U_sb_11);
    const ccomplex_t IT_0046 = IT_0021*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0041 + (-3)*IT_0044 + 3
      *IT_0047);
    const ccomplex_t IT_0049 = 0.166666666666667*IT_0048;
    const ccomplex_t IT_0050 = IT_0038*IT_0049;
    const ccomplex_t IT_0051 = IT_0036*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0011*IT_0013*IT_0015;
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = N_B4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0057 = IT_0026*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = N_W4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0060 = IT_0015*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = m_b*N_d4*e_em*IT_0011*conj(U_sb_11);
    const ccomplex_t IT_0063 = IT_0021*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0058 + (-3)*IT_0061 + 3
      *IT_0064);
    const ccomplex_t IT_0066 = 0.166666666666667*IT_0065;
    const ccomplex_t IT_0067 = IT_0055*IT_0066;
    const ccomplex_t IT_0068 = IT_0053*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0011*IT_0013*IT_0015;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0074 = IT_0026*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0077 = IT_0015*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = m_b*N_d3*e_em*IT_0011*conj(U_sb_11);
    const ccomplex_t IT_0080 = IT_0021*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0075 + (-3)*IT_0078 + 3
      *IT_0081);
    const ccomplex_t IT_0083 = 0.166666666666667*IT_0082;
    const ccomplex_t IT_0084 = IT_0072*IT_0083;
    const ccomplex_t IT_0085 = IT_0070*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = -IT_0035 + -IT_0052 + -IT_0069 + -IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0026;
    const ccomplex_t IT_0089 = (-0.333333333333333)*IT_0088;
    const ccomplex_t IT_0090 = conj(N_B1)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0091 = IT_0026*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = m_b*conj(N_d1)*e_em*IT_0011*conj(U_sb_01);
    const ccomplex_t IT_0094 = IT_0021*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0092 + 1.5*IT_0095);
    const ccomplex_t IT_0097 = (-0.333333333333333)*IT_0096;
    const ccomplex_t IT_0098 = IT_0089*IT_0097;
    const ccomplex_t IT_0099 = IT_0010*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0026;
    const ccomplex_t IT_0102 = (-0.333333333333333)*IT_0101;
    const ccomplex_t IT_0103 = conj(N_B4)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0104 = IT_0026*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = m_b*conj(N_d4)*e_em*IT_0011*conj(U_sb_01);
    const ccomplex_t IT_0107 = IT_0021*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + 1.5*IT_0108);
    const ccomplex_t IT_0110 = (-0.333333333333333)*IT_0109;
    const ccomplex_t IT_0111 = IT_0102*IT_0110;
    const ccomplex_t IT_0112 = IT_0053*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0026;
    const ccomplex_t IT_0115 = (-0.333333333333333)*IT_0114;
    const ccomplex_t IT_0116 = conj(N_B2)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0117 = IT_0026*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = m_b*conj(N_d2)*e_em*IT_0011*conj(U_sb_01);
    const ccomplex_t IT_0120 = IT_0021*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0118 + 1.5*IT_0121);
    const ccomplex_t IT_0123 = (-0.333333333333333)*IT_0122;
    const ccomplex_t IT_0124 = IT_0115*IT_0123;
    const ccomplex_t IT_0125 = IT_0036*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0026;
    const ccomplex_t IT_0128 = (-0.333333333333333)*IT_0127;
    const ccomplex_t IT_0129 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0130 = IT_0026*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = m_b*conj(N_d3)*e_em*IT_0011*conj(U_sb_01);
    const ccomplex_t IT_0133 = IT_0021*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0131 + 1.5*IT_0134);
    const ccomplex_t IT_0136 = (-0.333333333333333)*IT_0135;
    const ccomplex_t IT_0137 = IT_0128*IT_0136;
    const ccomplex_t IT_0138 = IT_0070*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = 0.166666666666667*IT_0008 + -IT_0100 + -IT_0113
       + -IT_0126 + -IT_0139;
    const ccomplex_t IT_0141 = s_34*IT_0002;
    const ccomplex_t IT_0142 = (-18)*IT_0141;
    const ccomplex_t IT_0143 = s_13*s_14;
    const ccomplex_t IT_0144 = 36*IT_0143;
    const ccomplex_t IT_0145 = IT_0142 + IT_0144;
    const ccomplex_t IT_0146 = m_s*s_14;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_01);
    const ccomplex_t IT_0148 = IT_0004*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0150 = IT_0003*IT_0148*IT_0149;
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0153 = IT_0032*IT_0089;
    const ccomplex_t IT_0154 = IT_0010*IT_0152*IT_0153;
    const ccomplex_t IT_0155 = IT_0066*IT_0102;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0157 = IT_0053*IT_0155*IT_0156;
    const ccomplex_t IT_0158 = m_s*IT_0035;
    const ccomplex_t IT_0159 = m_s*IT_0069;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0161 = IT_0049*IT_0115;
    const ccomplex_t IT_0162 = IT_0036*IT_0160*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0164 = IT_0083*IT_0128;
    const ccomplex_t IT_0165 = IT_0070*IT_0163*IT_0164;
    const ccomplex_t IT_0166 = m_s*IT_0052;
    const ccomplex_t IT_0167 = m_s*IT_0086;
    const ccomplex_t IT_0168 = 0.166666666666667*IT_0150 + -IT_0154 + -IT_0157
       + IT_0158 + IT_0159 + -IT_0162 + -IT_0165 + IT_0166 + IT_0167;
    const ccomplex_t IT_0169 = m_s*IT_0008;
    const ccomplex_t IT_0170 = 0.5*IT_0169;
    const ccomplex_t IT_0171 = m_s*IT_0113;
    const ccomplex_t IT_0172 = IT_0072*IT_0136;
    const ccomplex_t IT_0173 = IT_0070*IT_0163*IT_0172;
    const ccomplex_t IT_0174 = IT_0038*IT_0123;
    const ccomplex_t IT_0175 = IT_0036*IT_0160*IT_0174;
    const ccomplex_t IT_0176 = m_s*IT_0126;
    const ccomplex_t IT_0177 = m_s*IT_0139;
    const ccomplex_t IT_0178 = IT_0017*IT_0097;
    const ccomplex_t IT_0179 = IT_0010*IT_0152*IT_0178;
    const ccomplex_t IT_0180 = IT_0055*IT_0110;
    const ccomplex_t IT_0181 = IT_0053*IT_0156*IT_0180;
    const ccomplex_t IT_0182 = m_s*IT_0100;
    const ccomplex_t IT_0183 = (-0.166666666666667)*IT_0169 + IT_0171 + 
      -IT_0173 + -IT_0175 + IT_0176 + IT_0177 + -IT_0179 + -IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = 3*IT_0140;
    const ccomplex_t IT_0185 = IT_0009 + IT_0184;
    const ccomplex_t IT_0186 = 3*conj(IT_0140);
    const ccomplex_t IT_0187 = conj(IT_0009) + IT_0186;
    const ccomplex_t IT_0188 = conj(IT_0087)*(IT_0151 + 3*IT_0168) + IT_0087*
      (conj(IT_0151) + 3*conj(IT_0168)) + 0.166666666666667*(18*conj(IT_0009) +
       6*conj(IT_0140))*IT_0170 + 0.166666666666667*(18*IT_0009 + 6*IT_0140)
      *conj(IT_0170) + conj(IT_0183)*IT_0185 + IT_0183*IT_0187;
    const ccomplex_t IT_0189 = m_b*m_s*IT_0002;
    const ccomplex_t IT_0190 = 6*s_34;
    const ccomplex_t IT_0191 = 18*s_34;
    const ccomplex_t IT_0192 = m_b*m_s;
    const ccomplex_t IT_0193 = (-18)*conj(IT_0151);
    const ccomplex_t IT_0194 = (-6)*conj(IT_0168);
    const ccomplex_t IT_0195 = IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = (-18)*conj(IT_0168);
    const ccomplex_t IT_0197 = (-6)*conj(IT_0151);
    const ccomplex_t IT_0198 = IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = (-18)*IT_0168;
    const ccomplex_t IT_0200 = (-6)*IT_0168;
    const ccomplex_t IT_0201 = m_b*s_13;
    const ccomplex_t IT_0202 = (-0.166666666666667)*((-18)*conj(IT_0009) + (-6
      )*conj(IT_0140))*IT_0151 + conj(IT_0087)*(IT_0170 + 3*IT_0183) + IT_0087*
      (conj(IT_0170) + 3*conj(IT_0183)) + (-0.166666666666667)*IT_0009*IT_0195 +
       (-0.166666666666667)*IT_0140*IT_0198 + (-0.166666666666667)*conj(IT_0140)
      *IT_0199 + (-0.166666666666667)*conj(IT_0009)*IT_0200;
    const ccomplex_t IT_0203 = 12*IT_0143;
    const ccomplex_t IT_0204 = (-6)*IT_0141;
    const ccomplex_t IT_0205 = (IT_0009*conj(IT_0009) + IT_0087*conj(IT_0087) 
      + IT_0140*conj(IT_0140))*IT_0145 + 6*IT_0146*IT_0188 + (-6)*(conj(IT_0087)
      *IT_0185 + IT_0087*IT_0187)*IT_0189 + (conj(IT_0151)*IT_0168 + IT_0151
      *conj(IT_0168) + conj(IT_0170)*IT_0183 + IT_0170*conj(IT_0183))*IT_0190 + 
      (IT_0151*conj(IT_0151) + IT_0168*conj(IT_0168) + IT_0170*conj(IT_0170) +
       IT_0183*conj(IT_0183))*IT_0191 + IT_0192*(IT_0151*((-18)*conj(IT_0170) + 
      (-6)*conj(IT_0183)) + IT_0170*IT_0195 + IT_0183*IT_0198 + conj(IT_0183)
      *IT_0199 + conj(IT_0170)*IT_0200) + (-6)*IT_0201*IT_0202 + (conj(IT_0009)
      *IT_0140 + IT_0009*conj(IT_0140))*(IT_0203 + IT_0204);
    return create_ccomplex_return(IT_0205);
}

