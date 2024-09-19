#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_sb_1_to_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_sb_1_to_b_b(
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
    const creal_t m_sb_1 = param->m_sb_1;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0001 = cpow(IT_0000, 2);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*m_sG*IT_0001;
    const ccomplex_t IT_0003 = pow(m_b, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_sb_1, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = cpow((-2)*s_23 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = cos(beta);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = sin(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = pow(m_W, -1);
    const ccomplex_t IT_0020 = m_b*N_d1*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0013 + 1.5*IT_0022);
    const ccomplex_t IT_0024 = (-0.333333333333333)*IT_0023;
    const ccomplex_t IT_0025 = cpow(IT_0024, 2);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_N_1*IT_0025;
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0029 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0030 = IT_0010*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = m_b*N_d2*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0033 = IT_0018*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + 1.5*IT_0034);
    const ccomplex_t IT_0036 = (-0.333333333333333)*IT_0035;
    const ccomplex_t IT_0037 = cpow(IT_0036, 2);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_N_2*IT_0037;
    const ccomplex_t IT_0039 = IT_0028*IT_0038;
    const ccomplex_t IT_0040 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0041 = IT_0010*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = m_b*N_d4*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0044 = IT_0018*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 1.5*IT_0045);
    const ccomplex_t IT_0047 = (-0.333333333333333)*IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0047, 2);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*m_N_4*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_23 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = 0.166666666666667*IT_0007 + -IT_0027 + -IT_0039
       + -IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_13 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0055 = IT_0000*IT_0054;
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = m_b*IT_0057;
    const ccomplex_t IT_0059 = IT_0006*IT_0055;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = m_b*IT_0060;
    const ccomplex_t IT_0062 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0063 = IT_0010*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0066 = IT_0017*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = m_b*conj(N_d3)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0069 = IT_0018*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0064 + (-3)*IT_0067 + 3
      *IT_0070);
    const ccomplex_t IT_0072 = 0.166666666666667*IT_0071;
    const ccomplex_t IT_0073 = cpow(IT_0072, 2);
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*m_N_3*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_23 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0078 = IT_0010*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0081 = IT_0017*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = m_b*conj(N_d1)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0084 = IT_0018*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0079 + (-3)*IT_0082 + 3
      *IT_0085);
    const ccomplex_t IT_0087 = 0.166666666666667*IT_0086;
    const ccomplex_t IT_0088 = IT_0024*IT_0087;
    const ccomplex_t IT_0089 = IT_0008*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = m_b*IT_0090;
    const ccomplex_t IT_0092 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0093 = IT_0010*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0096 = IT_0017*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = m_b*conj(N_d2)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0099 = IT_0018*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0094 + (-3)*IT_0097 + 3
      *IT_0100);
    const ccomplex_t IT_0102 = 0.166666666666667*IT_0101;
    const ccomplex_t IT_0103 = IT_0036*IT_0102;
    const ccomplex_t IT_0104 = IT_0028*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = m_b*IT_0105;
    const ccomplex_t IT_0107 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0108 = IT_0010*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = m_b*N_d3*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0111 = IT_0018*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + 1.5*IT_0112);
    const ccomplex_t IT_0114 = (-0.333333333333333)*IT_0113;
    const ccomplex_t IT_0115 = IT_0072*IT_0114;
    const ccomplex_t IT_0116 = IT_0075*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = m_b*IT_0117;
    const ccomplex_t IT_0119 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0120 = IT_0010*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0123 = IT_0017*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = m_b*conj(N_d4)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0126 = IT_0018*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0121 + (-3)*IT_0124 + 3
      *IT_0127);
    const ccomplex_t IT_0129 = 0.166666666666667*IT_0128;
    const ccomplex_t IT_0130 = IT_0047*IT_0129;
    const ccomplex_t IT_0131 = IT_0050*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = m_b*IT_0132;
    const ccomplex_t IT_0134 = cpow(IT_0102, 2);
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*m_N_2*IT_0134;
    const ccomplex_t IT_0136 = IT_0028*IT_0135;
    const ccomplex_t IT_0137 = cpow(IT_0087, 2);
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*m_N_1*IT_0137;
    const ccomplex_t IT_0139 = IT_0008*IT_0138;
    const ccomplex_t IT_0140 = (-0.5)*IT_0058 + 0.166666666666667*IT_0061 +
       IT_0076 + -IT_0091 + -IT_0106 + -IT_0118 + -IT_0133 + IT_0136 + IT_0139;
    const ccomplex_t IT_0141 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0142 = IT_0026*IT_0141;
    const ccomplex_t IT_0143 = cpow(IT_0114, 2);
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*m_N_3*IT_0143;
    const ccomplex_t IT_0145 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0146 = IT_0144*IT_0145;
    const ccomplex_t IT_0147 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0148 = IT_0049*IT_0147;
    const ccomplex_t IT_0149 = (-0.5)*IT_0007 + -IT_0142 + -IT_0146 + -IT_0148;
    const ccomplex_t IT_0150 = 6*s_34;
    const ccomplex_t IT_0151 = IT_0130*IT_0147;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = IT_0088*IT_0141;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0156 = IT_0103*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = IT_0115*IT_0145;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = (-3)*IT_0060 + (-6)*IT_0152 + (-6)*IT_0154 + (
      -6)*IT_0157 + (-6)*IT_0159;
    const ccomplex_t IT_0161 = IT_0057 + IT_0160;
    const ccomplex_t IT_0162 = m_b*IT_0161;
    const ccomplex_t IT_0163 = 0.166666666666667*IT_0162;
    const ccomplex_t IT_0164 = 6*IT_0003;
    const ccomplex_t IT_0165 = 0.166666666666667*IT_0057 + 0.5*IT_0060 + 
      -IT_0152 + -IT_0154 + -IT_0157 + -IT_0159;
    const ccomplex_t IT_0166 = m_b*s_23;
    const ccomplex_t IT_0167 = 6*IT_0166;
    const ccomplex_t IT_0168 = (-0.5)*IT_0057 + (-0.166666666666667)*IT_0060 +
       IT_0090 + IT_0105 + IT_0117 + IT_0132;
    const ccomplex_t IT_0169 = 18*IT_0166;
    const ccomplex_t IT_0170 = IT_0002*IT_0053;
    const ccomplex_t IT_0171 = IT_0075*IT_0144;
    const ccomplex_t IT_0172 = (-0.5)*IT_0058 + 0.166666666666667*IT_0061 + 
      -IT_0091 + -IT_0106 + -IT_0118 + -IT_0133 + 0.5*IT_0170 + IT_0171;
    const ccomplex_t IT_0173 = (-18)*s_34;
    const ccomplex_t IT_0174 = cpow(IT_0054, 2);
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*m_sG*IT_0174;
    const ccomplex_t IT_0176 = IT_0006*IT_0175;
    const ccomplex_t IT_0177 = IT_0053*IT_0175;
    const ccomplex_t IT_0178 = cpow(IT_0129, 2);
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*m_N_4*IT_0178;
    const ccomplex_t IT_0180 = IT_0050*IT_0179;
    const ccomplex_t IT_0181 = 0.166666666666667*IT_0176 + (-0.5)*IT_0177 + 
      -IT_0180;
    const ccomplex_t IT_0182 = (-18)*IT_0003;
    const ccomplex_t IT_0183 = m_b*s_24;
    const ccomplex_t IT_0184 = (-18)*IT_0183;
    const ccomplex_t IT_0185 = IT_0038*IT_0155;
    const ccomplex_t IT_0186 = m_b*IT_0154;
    const ccomplex_t IT_0187 = m_b*IT_0152;
    const ccomplex_t IT_0188 = m_b*IT_0157;
    const ccomplex_t IT_0189 = m_b*IT_0159;
    const ccomplex_t IT_0190 = 0.166666666666667*IT_0058 + (-0.5)*IT_0061 + (
      -0.166666666666667)*IT_0170 + IT_0185 + -IT_0186 + -IT_0187 + -IT_0188 + 
      -IT_0189;
    const ccomplex_t IT_0191 = (-6)*s_34;
    const ccomplex_t IT_0192 = IT_0074*IT_0145;
    const ccomplex_t IT_0193 = IT_0147*IT_0179;
    const ccomplex_t IT_0194 = IT_0135*IT_0155;
    const ccomplex_t IT_0195 = IT_0138*IT_0141;
    const ccomplex_t IT_0196 = (-0.5)*IT_0176 + 0.166666666666667*IT_0177 + 
      -IT_0192 + -IT_0193 + -IT_0194 + -IT_0195;
    const ccomplex_t IT_0197 = (-6)*IT_0003;
    const ccomplex_t IT_0198 = (-6)*IT_0183;
    const ccomplex_t IT_0199 = 18*s_34*conj(IT_0052) + 18*IT_0003*conj(IT_0140
      ) + conj(IT_0149)*IT_0150 + conj(IT_0163)*IT_0164 + conj(IT_0165)*IT_0167 
      + conj(IT_0168)*IT_0169 + conj(IT_0172)*IT_0173 + conj(IT_0181)*IT_0182 +
       conj(IT_0168)*IT_0184 + conj(IT_0190)*IT_0191 + conj(IT_0196)*IT_0197 +
       conj(IT_0165)*IT_0198;
    const ccomplex_t IT_0200 = 18*s_34*conj(IT_0149) + conj(IT_0052)*IT_0150 +
       18*IT_0003*conj(IT_0163) + conj(IT_0140)*IT_0164 + IT_0167*conj(IT_0168) 
      + conj(IT_0165)*IT_0169 + conj(IT_0165)*IT_0184 + IT_0173*conj(IT_0190) +
       conj(IT_0172)*IT_0191 + IT_0182*conj(IT_0196) + conj(IT_0181)*IT_0197 +
       conj(IT_0168)*IT_0198;
    const ccomplex_t IT_0201 = conj(IT_0165)*IT_0167 + conj(IT_0168)*IT_0169 +
       18*IT_0003*conj(IT_0172) + conj(IT_0140)*IT_0173 + 18*s_34*conj(IT_0181) 
      + conj(IT_0052)*IT_0182 + conj(IT_0168)*IT_0184 + IT_0164*conj(IT_0190) +
       conj(IT_0163)*IT_0191 + IT_0150*conj(IT_0196) + conj(IT_0149)*IT_0197 +
       conj(IT_0165)*IT_0198;
    const ccomplex_t IT_0202 = IT_0167*conj(IT_0168) + conj(IT_0165)*IT_0169 +
       IT_0164*conj(IT_0172) + conj(IT_0163)*IT_0173 + IT_0150*conj(IT_0181) +
       conj(IT_0149)*IT_0182 + conj(IT_0165)*IT_0184 + 18*IT_0003*conj(IT_0190) 
      + conj(IT_0140)*IT_0191 + 18*s_34*conj(IT_0196) + conj(IT_0052)*IT_0197 +
       conj(IT_0168)*IT_0198;
    const ccomplex_t IT_0203 = 6*IT_0183;
    const ccomplex_t IT_0204 = 18*IT_0183;
    const ccomplex_t IT_0205 = (-18)*IT_0166;
    const ccomplex_t IT_0206 = (-6)*IT_0166;
    const ccomplex_t IT_0207 = 18*IT_0003*conj(IT_0052) + 18*s_34*conj(IT_0140
      ) + IT_0150*conj(IT_0163) + conj(IT_0149)*IT_0164 + IT_0173*conj(IT_0181) 
      + conj(IT_0172)*IT_0182 + IT_0191*conj(IT_0196) + conj(IT_0190)*IT_0197 +
       conj(IT_0165)*IT_0203 + conj(IT_0168)*IT_0204 + conj(IT_0168)*IT_0205 +
       conj(IT_0165)*IT_0206;
    const ccomplex_t IT_0208 = 18*IT_0003*conj(IT_0149) + conj(IT_0140)
      *IT_0150 + 18*s_34*conj(IT_0163) + conj(IT_0052)*IT_0164 + IT_0182*conj
      (IT_0190) + conj(IT_0181)*IT_0191 + IT_0173*conj(IT_0196) + conj(IT_0172)
      *IT_0197 + conj(IT_0168)*IT_0203 + conj(IT_0165)*IT_0204 + conj(IT_0165)
      *IT_0205 + conj(IT_0168)*IT_0206;
    const ccomplex_t IT_0209 = 18*s_34*conj(IT_0172) + conj(IT_0052)*IT_0173 +
       18*IT_0003*conj(IT_0181) + conj(IT_0140)*IT_0182 + IT_0150*conj(IT_0190) 
      + conj(IT_0149)*IT_0191 + IT_0164*conj(IT_0196) + conj(IT_0163)*IT_0197 +
       conj(IT_0165)*IT_0203 + conj(IT_0168)*IT_0204 + conj(IT_0168)*IT_0205 +
       conj(IT_0165)*IT_0206;
    const ccomplex_t IT_0210 = IT_0150*conj(IT_0172) + conj(IT_0149)*IT_0173 +
       IT_0164*conj(IT_0181) + conj(IT_0163)*IT_0182 + 18*s_34*conj(IT_0190) +
       conj(IT_0052)*IT_0191 + 18*IT_0003*conj(IT_0196) + conj(IT_0140)*IT_0197 
      + conj(IT_0168)*IT_0203 + conj(IT_0165)*IT_0204 + conj(IT_0165)*IT_0205 +
       conj(IT_0168)*IT_0206;
    const ccomplex_t IT_0211 = s_34*IT_0005;
    const ccomplex_t IT_0212 = (-6)*IT_0211;
    const ccomplex_t IT_0213 = s_23*s_24;
    const ccomplex_t IT_0214 = 12*IT_0213;
    const ccomplex_t IT_0215 = IT_0212 + IT_0214;
    const ccomplex_t IT_0216 = (-18)*IT_0211;
    const ccomplex_t IT_0217 = 36*IT_0213;
    const ccomplex_t IT_0218 = IT_0216 + IT_0217;
    const ccomplex_t IT_0219 = IT_0003*IT_0005;
    const ccomplex_t IT_0220 = (-18)*IT_0219;
    const ccomplex_t IT_0221 = (-6)*IT_0219;
    const ccomplex_t IT_0222 = conj(IT_0149)*IT_0167 + conj(IT_0052)*IT_0169 +
       conj(IT_0181)*IT_0184 + conj(IT_0196)*IT_0198 + conj(IT_0163)*IT_0203 +
       conj(IT_0140)*IT_0204 + conj(IT_0172)*IT_0205 + conj(IT_0190)*IT_0206 +
       conj(IT_0165)*IT_0215 + conj(IT_0168)*IT_0218 + conj(IT_0168)*IT_0220 +
       conj(IT_0165)*IT_0221;
    const ccomplex_t IT_0223 = conj(IT_0052)*IT_0167 + conj(IT_0149)*IT_0169 +
       IT_0184*conj(IT_0196) + conj(IT_0181)*IT_0198 + conj(IT_0140)*IT_0203 +
       conj(IT_0163)*IT_0204 + conj(IT_0190)*IT_0205 + conj(IT_0172)*IT_0206 +
       conj(IT_0168)*IT_0215 + conj(IT_0165)*IT_0218 + conj(IT_0165)*IT_0220 +
       conj(IT_0168)*IT_0221;
    const ccomplex_t IT_0224 = IT_0169*conj(IT_0181) + conj(IT_0052)*IT_0184 +
       IT_0167*conj(IT_0196) + conj(IT_0149)*IT_0198 + conj(IT_0190)*IT_0203 +
       conj(IT_0172)*IT_0204 + conj(IT_0140)*IT_0205 + conj(IT_0163)*IT_0206 +
       conj(IT_0165)*IT_0215 + conj(IT_0168)*IT_0218 + conj(IT_0168)*IT_0220 +
       conj(IT_0165)*IT_0221;
    const ccomplex_t IT_0225 = IT_0167*conj(IT_0181) + conj(IT_0149)*IT_0184 +
       IT_0169*conj(IT_0196) + conj(IT_0052)*IT_0198 + conj(IT_0172)*IT_0203 +
       conj(IT_0190)*IT_0204 + conj(IT_0163)*IT_0205 + conj(IT_0140)*IT_0206 +
       conj(IT_0168)*IT_0215 + conj(IT_0165)*IT_0218 + conj(IT_0165)*IT_0220 +
       conj(IT_0168)*IT_0221;
    const ccomplex_t IT_0226 = IT_0052*IT_0199 + IT_0149*IT_0200 + IT_0181
      *IT_0201 + IT_0196*IT_0202 + IT_0140*IT_0207 + IT_0163*IT_0208 + IT_0172
      *IT_0209 + IT_0190*IT_0210 + IT_0168*IT_0222 + IT_0165*IT_0223 + IT_0168
      *IT_0224 + IT_0165*IT_0225;
    return create_ccomplex_return(IT_0226);
}

