#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
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
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cpow(V_tb, 2);
    const ccomplex_t IT_0009 = m_b*mu_h*IT_0000;
    const ccomplex_t IT_0010 = sin(alpha);
    const ccomplex_t IT_0011 = sin(beta);
    const ccomplex_t IT_0012 = -IT_0010*IT_0011;
    const ccomplex_t IT_0013 = IT_0001*IT_0003;
    const ccomplex_t IT_0014 = IT_0012 + IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -2);
    const ccomplex_t IT_0017 = pow(m_b, 2);
    const ccomplex_t IT_0018 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0014*(IT_0005 + 1.
      /3*IT_0004*IT_0016) + IT_0000*IT_0002*IT_0003*IT_0005*IT_0017)*conj
      (U_sb_01) + -1./2*IT_0002*IT_0005*IT_0009*IT_0010*conj(U_sb_11)) + U_sb_10
      *(IT_0002*IT_0005*IT_0009*IT_0010*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0004*IT_0014*IT_0016 + (-3)*IT_0000*IT_0002*IT_0003*IT_0005*IT_0017)
      *conj(U_sb_11));
    const ccomplex_t IT_0019 = (0 + _Complex_I*-1)*(e_em*IT_0018 + (
      -1.4142135623731)*A_b*IT_0003*IT_0008*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = pow(m_sb_1, 2);
    const ccomplex_t IT_0023 = pow(m_sb_2, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + -IT_0022 + -IT_0023 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0005*IT_0010;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0001*IT_0010;
    const ccomplex_t IT_0030 = IT_0003*IT_0011;
    const ccomplex_t IT_0031 = IT_0029 + IT_0030;
    const ccomplex_t IT_0032 = 2*U_sb_00*((IT_0000*IT_0002*IT_0005*IT_0010
      *IT_0017 + (-0.5)*m_W*(IT_0005 + 1./3*IT_0004*IT_0016)*IT_0031)*conj
      (U_sb_01) + 1./2*IT_0002*IT_0003*IT_0005*IT_0009*conj(U_sb_11)) + U_sb_10*
      (IT_0002*IT_0003*IT_0005*IT_0009*conj(U_sb_01) + 2*(IT_0000*IT_0002
      *IT_0005*IT_0010*IT_0017 + -1./3*m_W*IT_0004*IT_0016*IT_0031)*conj(U_sb_11
      ));
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(e_em*IT_0032 +
       1.4142135623731*A_b*IT_0008*IT_0010*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0028*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0022 + -IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = -IT_0026 + -IT_0038;
    const ccomplex_t IT_0040 = 36*s_34;
    const ccomplex_t IT_0041 = cpow(IT_0015, -1);
    const ccomplex_t IT_0042 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = IT_0002*IT_0005;
    const ccomplex_t IT_0046 = m_b*N_d1*e_em*IT_0000*U_sb_00;
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0044 + 1.5*IT_0048);
    const ccomplex_t IT_0050 = (-0.333333333333333)*IT_0049;
    const ccomplex_t IT_0051 = N_B1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0052 = IT_0041*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = N_W1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0055 = IT_0005*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = m_b*N_d1*e_em*IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0058 = IT_0045*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0053 + (-3)*IT_0056 + 3
      *IT_0059);
    const ccomplex_t IT_0061 = 0.166666666666667*IT_0060;
    const ccomplex_t IT_0062 = IT_0050*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0064 = cpow((-2)*s_23 + IT_0017 + IT_0023 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0062*IT_0063*IT_0064;
    const ccomplex_t IT_0066 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0067 = IT_0041*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = m_b*N_d2*e_em*IT_0000*U_sb_00;
    const ccomplex_t IT_0070 = IT_0045*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0068 + 1.5*IT_0071);
    const ccomplex_t IT_0073 = (-0.333333333333333)*IT_0072;
    const ccomplex_t IT_0074 = N_B2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0075 = IT_0041*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = N_W2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0078 = IT_0005*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = m_b*N_d2*e_em*IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0081 = IT_0045*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0076 + (-3)*IT_0079 + 3
      *IT_0082);
    const ccomplex_t IT_0084 = 0.166666666666667*IT_0083;
    const ccomplex_t IT_0085 = IT_0073*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0087 = cpow((-2)*s_23 + IT_0017 + IT_0023 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0085*IT_0086*IT_0087;
    const ccomplex_t IT_0089 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0090 = IT_0041*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = m_b*N_d3*e_em*IT_0000*U_sb_00;
    const ccomplex_t IT_0093 = IT_0045*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0091 + 1.5*IT_0094);
    const ccomplex_t IT_0096 = (-0.333333333333333)*IT_0095;
    const ccomplex_t IT_0097 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0098 = IT_0041*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0101 = IT_0005*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = m_b*N_d3*e_em*IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0104 = IT_0045*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(IT_0099 + (-3)*IT_0102 + 3
      *IT_0105);
    const ccomplex_t IT_0107 = 0.166666666666667*IT_0106;
    const ccomplex_t IT_0108 = IT_0096*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0110 = cpow((-2)*s_23 + IT_0017 + IT_0023 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0108*IT_0109*IT_0110;
    const ccomplex_t IT_0112 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0113 = IT_0041*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = m_b*N_d4*e_em*IT_0000*U_sb_00;
    const ccomplex_t IT_0116 = IT_0045*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0114 + 1.5*IT_0117);
    const ccomplex_t IT_0119 = (-0.333333333333333)*IT_0118;
    const ccomplex_t IT_0120 = N_B4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0121 = IT_0041*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = N_W4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0124 = IT_0005*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = m_b*N_d4*e_em*IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0127 = IT_0045*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0122 + (-3)*IT_0125 + 3
      *IT_0128);
    const ccomplex_t IT_0130 = 0.166666666666667*IT_0129;
    const ccomplex_t IT_0131 = IT_0119*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0133 = cpow((-2)*s_23 + IT_0017 + IT_0023 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0134 = IT_0131*IT_0132*IT_0133;
    const ccomplex_t IT_0135 = pow(m_sG, 2);
    const ccomplex_t IT_0136 = cpow((-2)*s_23 + IT_0017 + IT_0023 + -IT_0135 +
       reg_prop, -1);
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_01);
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = IT_0136*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = m_b*IT_0141;
    const ccomplex_t IT_0143 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0144 = IT_0041*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0147 = IT_0005*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = m_b*conj(N_d1)*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0150 = IT_0045*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0145 + (-3)*IT_0148 + 3
      *IT_0151);
    const ccomplex_t IT_0153 = 0.166666666666667*IT_0152;
    const ccomplex_t IT_0154 = IT_0061*IT_0153;
    const ccomplex_t IT_0155 = IT_0064*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = m_b*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0159 = IT_0041*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0162 = IT_0005*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = m_b*conj(N_d2)*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0165 = IT_0045*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0160 + (-3)*IT_0163 + 3
      *IT_0166);
    const ccomplex_t IT_0168 = 0.166666666666667*IT_0167;
    const ccomplex_t IT_0169 = IT_0084*IT_0168;
    const ccomplex_t IT_0170 = IT_0087*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = m_b*IT_0171;
    const ccomplex_t IT_0173 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0174 = IT_0041*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0177 = IT_0005*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = m_b*conj(N_d3)*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0180 = IT_0045*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*(IT_0175 + (-3)*IT_0178 + 3
      *IT_0181);
    const ccomplex_t IT_0183 = 0.166666666666667*IT_0182;
    const ccomplex_t IT_0184 = IT_0107*IT_0183;
    const ccomplex_t IT_0185 = IT_0110*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = m_b*IT_0186;
    const ccomplex_t IT_0188 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0189 = IT_0041*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0192 = IT_0005*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = m_b*conj(N_d4)*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0195 = IT_0045*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0190 + (-3)*IT_0193 + 3
      *IT_0196);
    const ccomplex_t IT_0198 = 0.166666666666667*IT_0197;
    const ccomplex_t IT_0199 = IT_0130*IT_0198;
    const ccomplex_t IT_0200 = IT_0133*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = m_b*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0204 = IT_0138*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0206 = IT_0136*IT_0204*IT_0205;
    const ccomplex_t IT_0207 = -IT_0065 + -IT_0088 + -IT_0111 + -IT_0134 + (
      -0.166666666666667)*IT_0142 + IT_0157 + IT_0172 + IT_0187 + IT_0202 +
       0.166666666666667*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_11);
    const ccomplex_t IT_0209 = IT_0137*IT_0208;
    const ccomplex_t IT_0210 = IT_0136*IT_0205*IT_0209;
    const ccomplex_t IT_0211 = conj(N_B1)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0212 = IT_0041*IT_0211;
    const ccomplex_t IT_0213 = 1.4142135623731*IT_0212;
    const ccomplex_t IT_0214 = m_b*conj(N_d1)*e_em*IT_0000*conj(U_sb_01);
    const ccomplex_t IT_0215 = IT_0045*IT_0214;
    const ccomplex_t IT_0216 = 1.4142135623731*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*(IT_0213 + 1.5*IT_0216);
    const ccomplex_t IT_0218 = (-0.333333333333333)*IT_0217;
    const ccomplex_t IT_0219 = IT_0153*IT_0218;
    const ccomplex_t IT_0220 = IT_0063*IT_0064*IT_0219;
    const ccomplex_t IT_0221 = conj(N_B2)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0222 = IT_0041*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = m_b*conj(N_d2)*e_em*IT_0000*conj(U_sb_01);
    const ccomplex_t IT_0225 = IT_0045*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0223 + 1.5*IT_0226);
    const ccomplex_t IT_0228 = (-0.333333333333333)*IT_0227;
    const ccomplex_t IT_0229 = IT_0168*IT_0228;
    const ccomplex_t IT_0230 = IT_0086*IT_0087*IT_0229;
    const ccomplex_t IT_0231 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0232 = IT_0041*IT_0231;
    const ccomplex_t IT_0233 = 1.4142135623731*IT_0232;
    const ccomplex_t IT_0234 = m_b*conj(N_d3)*e_em*IT_0000*conj(U_sb_01);
    const ccomplex_t IT_0235 = IT_0045*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*(IT_0233 + 1.5*IT_0236);
    const ccomplex_t IT_0238 = (-0.333333333333333)*IT_0237;
    const ccomplex_t IT_0239 = IT_0183*IT_0238;
    const ccomplex_t IT_0240 = IT_0109*IT_0110*IT_0239;
    const ccomplex_t IT_0241 = conj(N_B4)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0242 = IT_0041*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = m_b*conj(N_d4)*e_em*IT_0000*conj(U_sb_01);
    const ccomplex_t IT_0245 = IT_0045*IT_0244;
    const ccomplex_t IT_0246 = 1.4142135623731*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*(IT_0243 + 1.5*IT_0246);
    const ccomplex_t IT_0248 = (-0.333333333333333)*IT_0247;
    const ccomplex_t IT_0249 = IT_0198*IT_0248;
    const ccomplex_t IT_0250 = IT_0132*IT_0133*IT_0249;
    const ccomplex_t IT_0251 = IT_0119*IT_0248;
    const ccomplex_t IT_0252 = IT_0133*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = m_b*IT_0253;
    const ccomplex_t IT_0255 = IT_0073*IT_0228;
    const ccomplex_t IT_0256 = IT_0087*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = m_b*IT_0257;
    const ccomplex_t IT_0259 = IT_0096*IT_0238;
    const ccomplex_t IT_0260 = IT_0110*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*IT_0260;
    const ccomplex_t IT_0262 = m_b*IT_0261;
    const ccomplex_t IT_0263 = IT_0203*IT_0208;
    const ccomplex_t IT_0264 = IT_0136*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = m_b*IT_0265;
    const ccomplex_t IT_0267 = IT_0050*IT_0218;
    const ccomplex_t IT_0268 = IT_0064*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*IT_0268;
    const ccomplex_t IT_0270 = m_b*IT_0269;
    const ccomplex_t IT_0271 = 0.166666666666667*IT_0210 + -IT_0220 + -IT_0230
       + -IT_0240 + -IT_0250 + IT_0254 + IT_0258 + IT_0262 + (-0.166666666666667
      )*IT_0266 + IT_0270;
    const ccomplex_t IT_0272 = -IT_0017;
    const ccomplex_t IT_0273 = s_34 + IT_0272;
    const ccomplex_t IT_0274 = 6*IT_0273;
    const ccomplex_t IT_0275 = 0.166666666666667*IT_0141 + -IT_0156 + -IT_0171
       + -IT_0186 + -IT_0201;
    const ccomplex_t IT_0276 = -IT_0253 + -IT_0257 + -IT_0261 +
       0.166666666666667*IT_0265 + -IT_0269;
    const ccomplex_t IT_0277 = -s_14;
    const ccomplex_t IT_0278 = s_13 + IT_0277;
    const ccomplex_t IT_0279 = m_b*IT_0278;
    const ccomplex_t IT_0280 = 6*IT_0279;
    const ccomplex_t IT_0281 = conj(IT_0039)*IT_0280;
    const ccomplex_t IT_0282 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0022 + -IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0283 = IT_0004*IT_0041;
    const ccomplex_t IT_0284 = e_em*IT_0283;
    const ccomplex_t IT_0285 = IT_0005*IT_0015;
    const ccomplex_t IT_0286 = e_em*IT_0285;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*(IT_0284 + 3*IT_0286);
    const ccomplex_t IT_0288 = (-0.166666666666667)*IT_0287;
    const ccomplex_t IT_0289 = IT_0282*IT_0288;
    const ccomplex_t IT_0290 = (0 + _Complex_I*1)*IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0005
      *IT_0015 + 1./3*IT_0004*IT_0041)*U_sb_00*conj(U_sb_01) + IT_0004*IT_0041
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0292 = (-0.666666666666667)*IT_0291;
    const ccomplex_t IT_0293 = IT_0290*IT_0292;
    const ccomplex_t IT_0294 = (-0.5)*IT_0141 + -IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0041;
    const ccomplex_t IT_0296 = 0.333333333333333*IT_0295;
    const ccomplex_t IT_0297 = IT_0282*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = IT_0292*IT_0298;
    const ccomplex_t IT_0300 = (-0.5)*IT_0265 + -IT_0299;
    const ccomplex_t IT_0301 = 18*IT_0279;
    const ccomplex_t IT_0302 = conj(IT_0039)*IT_0301;
    const ccomplex_t IT_0303 = (-0.333333333333333)*IT_0291;
    const ccomplex_t IT_0304 = m_b*IT_0303;
    const ccomplex_t IT_0305 = m_b*IT_0292;
    const ccomplex_t IT_0306 = -IT_0305;
    const ccomplex_t IT_0307 = IT_0304 + IT_0306;
    const ccomplex_t IT_0308 = IT_0298*IT_0307;
    const ccomplex_t IT_0309 = IT_0290*IT_0307;
    const ccomplex_t IT_0310 = (-0.5)*IT_0210 + 0.5*IT_0266 + -IT_0308 +
       IT_0309;
    const ccomplex_t IT_0311 = (-18)*IT_0017;
    const ccomplex_t IT_0312 = 18*s_34;
    const ccomplex_t IT_0313 = IT_0311 + IT_0312;
    const ccomplex_t IT_0314 = m_b*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0315 = 0.5*IT_0314;
    const ccomplex_t IT_0316 = m_b*e_em*mu_h*IT_0000*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0317 = IT_0002*IT_0005*IT_0011;
    const ccomplex_t IT_0318 = IT_0316*IT_0317;
    const ccomplex_t IT_0319 = m_b*e_em*mu_h*IT_0000*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0320 = IT_0317*IT_0319;
    const ccomplex_t IT_0321 = A_b*IT_0008*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0322 = IT_0001*IT_0321;
    const ccomplex_t IT_0323 = 1.4142135623731*IT_0322;
    const ccomplex_t IT_0324 = A_b*IT_0008*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0325 = IT_0001*IT_0324;
    const ccomplex_t IT_0326 = 1.4142135623731*IT_0325;
    const ccomplex_t IT_0327 = (-0.5)*IT_0318 + 0.5*IT_0320 + 0.5*IT_0323 + (
      -0.5)*IT_0326;
    const ccomplex_t IT_0328 = IT_0315*IT_0327;
    const ccomplex_t IT_0329 = pow(m_Z, 2);
    const ccomplex_t IT_0330 = cpow((-2)*s_12 + -IT_0022 + -IT_0023 + IT_0329 
      + -reg_prop, -1);
    const ccomplex_t IT_0331 = IT_0328*IT_0330;
    const ccomplex_t IT_0332 = (0 + _Complex_I*1)*IT_0331;
    const ccomplex_t IT_0333 = m_b*e_em*IT_0000*IT_0002*IT_0005*IT_0011;
    const ccomplex_t IT_0334 = (-0.5)*IT_0333;
    const ccomplex_t IT_0335 = IT_0005*IT_0316;
    const ccomplex_t IT_0336 = IT_0005*IT_0319;
    const ccomplex_t IT_0337 = IT_0011*IT_0321;
    const ccomplex_t IT_0338 = 1.4142135623731*IT_0337;
    const ccomplex_t IT_0339 = IT_0011*IT_0324;
    const ccomplex_t IT_0340 = 1.4142135623731*IT_0339;
    const ccomplex_t IT_0341 = (-0.5)*IT_0335 + 0.5*IT_0336 + (-0.5)*IT_0338 +
       0.5*IT_0340;
    const ccomplex_t IT_0342 = IT_0334*IT_0341;
    const ccomplex_t IT_0343 = cpow((-2)*s_12 + -IT_0022 + -IT_0023 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0344 = IT_0342*IT_0343;
    const ccomplex_t IT_0345 = (0 + _Complex_I*1)*IT_0344;
    const ccomplex_t IT_0346 = -IT_0332 + -IT_0345;
    const ccomplex_t IT_0347 = 18*IT_0017;
    const ccomplex_t IT_0348 = IT_0312 + IT_0347;
    const ccomplex_t IT_0349 = 0.5*IT_0142 + (-0.5)*IT_0206 + IT_0308 + 
      -IT_0309;
    const ccomplex_t IT_0350 = 6*s_34;
    const ccomplex_t IT_0351 = m_b*s_13;
    const ccomplex_t IT_0352 = 6*IT_0351;
    const ccomplex_t IT_0353 = 18*IT_0351;
    const ccomplex_t IT_0354 = m_b*s_14;
    const ccomplex_t IT_0355 = (-18)*IT_0354;
    const ccomplex_t IT_0356 = (-6)*IT_0017;
    const ccomplex_t IT_0357 = (-6)*IT_0354;
    const ccomplex_t IT_0358 = conj(IT_0310)*(IT_0310*IT_0312 + IT_0039
      *IT_0313 + IT_0346*IT_0348 + IT_0311*IT_0349 + IT_0271*IT_0350 + IT_0276
      *IT_0352 + IT_0300*IT_0353 + IT_0294*IT_0355 + IT_0207*IT_0356 + IT_0275
      *IT_0357);
    const ccomplex_t IT_0359 = (-18)*s_34 + IT_0311;
    const ccomplex_t IT_0360 = conj(IT_0349)*(IT_0310*IT_0311 + IT_0039
      *IT_0313 + IT_0312*IT_0349 + IT_0207*IT_0350 + IT_0275*IT_0352 + IT_0294
      *IT_0353 + IT_0300*IT_0355 + IT_0271*IT_0356 + IT_0276*IT_0357 + IT_0346
      *IT_0359);
    const ccomplex_t IT_0361 = s_34 + IT_0017;
    const ccomplex_t IT_0362 = 6*IT_0361;
    const ccomplex_t IT_0363 = s_13 + s_14;
    const ccomplex_t IT_0364 = m_b*IT_0363;
    const ccomplex_t IT_0365 = 6*IT_0364;
    const ccomplex_t IT_0366 = 18*IT_0364;
    const ccomplex_t IT_0367 = (-18)*IT_0364;
    const ccomplex_t IT_0368 = (-6)*IT_0361;
    const ccomplex_t IT_0369 = (-6)*IT_0364;
    const ccomplex_t IT_0370 = conj(IT_0346)*((36*IT_0017 + IT_0040)*IT_0346 +
       IT_0310*IT_0348 + IT_0349*IT_0359 + IT_0271*IT_0362 + IT_0276*IT_0365 +
       IT_0300*IT_0366 + IT_0294*IT_0367 + IT_0207*IT_0368 + IT_0275*IT_0369);
    const ccomplex_t IT_0371 = conj(IT_0039)*IT_0313;
    const ccomplex_t IT_0372 = IT_0039*IT_0274;
    const ccomplex_t IT_0373 = s_34*IT_0022;
    const ccomplex_t IT_0374 = (-6)*IT_0373;
    const ccomplex_t IT_0375 = s_13*s_14;
    const ccomplex_t IT_0376 = 12*IT_0375;
    const ccomplex_t IT_0377 = IT_0374 + IT_0376;
    const ccomplex_t IT_0378 = (-18)*IT_0373;
    const ccomplex_t IT_0379 = 36*IT_0375;
    const ccomplex_t IT_0380 = IT_0378 + IT_0379;
    const ccomplex_t IT_0381 = IT_0017*IT_0022;
    const ccomplex_t IT_0382 = (-18)*IT_0381;
    const ccomplex_t IT_0383 = (-6)*IT_0381;
    const ccomplex_t IT_0384 = IT_0039*IT_0280;
    const ccomplex_t IT_0385 = IT_0039*IT_0301;
    const ccomplex_t IT_0386 = conj(IT_0039)*(IT_0039*((-36)*IT_0017 + IT_0040
      ) + (IT_0207 + IT_0271)*IT_0274) + (IT_0275 + IT_0276)*IT_0281 + (IT_0294 
      + IT_0300)*IT_0302 + IT_0358 + IT_0360 + IT_0370 + (IT_0310 + IT_0349)
      *IT_0371 + conj(IT_0271)*(IT_0207*IT_0311 + IT_0271*IT_0312 + IT_0310
      *IT_0350 + IT_0300*IT_0352 + IT_0276*IT_0353 + IT_0275*IT_0355 + IT_0349
      *IT_0356 + IT_0294*IT_0357 + IT_0346*IT_0362 + IT_0372) + conj(IT_0207)*
      (IT_0271*IT_0311 + IT_0207*IT_0312 + IT_0349*IT_0350 + IT_0294*IT_0352 +
       IT_0275*IT_0353 + IT_0276*IT_0355 + IT_0310*IT_0356 + IT_0300*IT_0357 +
       IT_0346*IT_0368 + IT_0372) + conj(IT_0276)*(IT_0310*IT_0352 + IT_0271
      *IT_0353 + IT_0207*IT_0355 + IT_0349*IT_0357 + IT_0346*IT_0365 + IT_0300
      *IT_0377 + IT_0276*IT_0380 + IT_0275*IT_0382 + IT_0294*IT_0383 + IT_0384) 
      + conj(IT_0275)*(IT_0349*IT_0352 + IT_0207*IT_0353 + IT_0271*IT_0355 +
       IT_0310*IT_0357 + IT_0346*IT_0369 + IT_0294*IT_0377 + IT_0275*IT_0380 +
       IT_0276*IT_0382 + IT_0300*IT_0383 + IT_0384) + conj(IT_0300)*(IT_0271
      *IT_0352 + IT_0310*IT_0353 + IT_0349*IT_0355 + IT_0207*IT_0357 + IT_0346
      *IT_0366 + IT_0276*IT_0377 + IT_0300*IT_0380 + IT_0294*IT_0382 + IT_0275
      *IT_0383 + IT_0385) + conj(IT_0294)*(IT_0207*IT_0352 + IT_0349*IT_0353 +
       IT_0310*IT_0355 + IT_0271*IT_0357 + IT_0346*IT_0367 + IT_0275*IT_0377 +
       IT_0294*IT_0380 + IT_0300*IT_0382 + IT_0276*IT_0383 + IT_0385);
    return create_ccomplex_return(IT_0386);
}

