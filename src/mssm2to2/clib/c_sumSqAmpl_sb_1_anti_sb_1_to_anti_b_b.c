#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_anti_b_b(
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
    const ccomplex_t U_sb_10 = param->U_sb_10;
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
      (U_sb_00) + -1./2*IT_0002*IT_0005*IT_0009*IT_0010*conj(U_sb_10)) + U_sb_10
      *(IT_0002*IT_0005*IT_0009*IT_0010*conj(U_sb_00) + 0.666666666666667*(m_W
      *IT_0004*IT_0014*IT_0016 + (-3)*IT_0000*IT_0002*IT_0003*IT_0005*IT_0017)
      *conj(U_sb_10));
    const ccomplex_t IT_0019 = (0 + _Complex_I*-1)*(e_em*IT_0018 + (
      -1.4142135623731)*A_b*IT_0003*IT_0008*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = pow(m_sb_1, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + (-2)*IT_0022 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0005*IT_0010;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0001*IT_0010;
    const ccomplex_t IT_0029 = IT_0003*IT_0011;
    const ccomplex_t IT_0030 = IT_0028 + IT_0029;
    const ccomplex_t IT_0031 = 2*U_sb_00*((IT_0000*IT_0002*IT_0005*IT_0010
      *IT_0017 + (-0.5)*m_W*(IT_0005 + 1./3*IT_0004*IT_0016)*IT_0030)*conj
      (U_sb_00) + 1./2*IT_0002*IT_0003*IT_0005*IT_0009*conj(U_sb_10)) + U_sb_10*
      (IT_0002*IT_0003*IT_0005*IT_0009*conj(U_sb_00) + 2*(IT_0000*IT_0002
      *IT_0005*IT_0010*IT_0017 + -1./3*m_W*IT_0004*IT_0016*IT_0030)*conj(U_sb_10
      ));
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(e_em*IT_0031 +
       1.4142135623731*A_b*IT_0008*IT_0010*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -IT_0025 + -IT_0037;
    const ccomplex_t IT_0039 = 36*s_34;
    const ccomplex_t IT_0040 = cpow(IT_0015, -1);
    const ccomplex_t IT_0041 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0045 = IT_0005*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = IT_0002*IT_0005;
    const ccomplex_t IT_0048 = m_b*N_d1*e_em*IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0043 + (-3)*IT_0046 + 3
      *IT_0050);
    const ccomplex_t IT_0052 = 0.166666666666667*IT_0051;
    const ccomplex_t IT_0053 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0054 = IT_0040*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = m_b*N_d1*e_em*IT_0000*U_sb_00;
    const ccomplex_t IT_0057 = IT_0047*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + 1.5*IT_0058);
    const ccomplex_t IT_0060 = (-0.333333333333333)*IT_0059;
    const ccomplex_t IT_0061 = IT_0052*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0063 = cpow((-2)*s_23 + IT_0017 + IT_0022 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0061*IT_0062*IT_0063;
    const ccomplex_t IT_0065 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0066 = IT_0040*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0069 = IT_0005*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = m_b*N_d2*e_em*IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0072 = IT_0047*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0067 + (-3)*IT_0070 + 3
      *IT_0073);
    const ccomplex_t IT_0075 = 0.166666666666667*IT_0074;
    const ccomplex_t IT_0076 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0077 = IT_0040*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = m_b*N_d2*e_em*IT_0000*U_sb_00;
    const ccomplex_t IT_0080 = IT_0047*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0078 + 1.5*IT_0081);
    const ccomplex_t IT_0083 = (-0.333333333333333)*IT_0082;
    const ccomplex_t IT_0084 = IT_0075*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0017 + IT_0022 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0084*IT_0085*IT_0086;
    const ccomplex_t IT_0088 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0089 = IT_0040*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0092 = IT_0005*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = m_b*N_d3*e_em*IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0095 = IT_0047*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0090 + (-3)*IT_0093 + 3
      *IT_0096);
    const ccomplex_t IT_0098 = 0.166666666666667*IT_0097;
    const ccomplex_t IT_0099 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0100 = IT_0040*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = m_b*N_d3*e_em*IT_0000*U_sb_00;
    const ccomplex_t IT_0103 = IT_0047*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + 1.5*IT_0104);
    const ccomplex_t IT_0106 = (-0.333333333333333)*IT_0105;
    const ccomplex_t IT_0107 = IT_0098*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0109 = cpow((-2)*s_23 + IT_0017 + IT_0022 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0110 = IT_0107*IT_0108*IT_0109;
    const ccomplex_t IT_0111 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0112 = IT_0040*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0115 = IT_0005*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = m_b*N_d4*e_em*IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0118 = IT_0047*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0113 + (-3)*IT_0116 + 3
      *IT_0119);
    const ccomplex_t IT_0121 = 0.166666666666667*IT_0120;
    const ccomplex_t IT_0122 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0123 = IT_0040*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = m_b*N_d4*e_em*IT_0000*U_sb_00;
    const ccomplex_t IT_0126 = IT_0047*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0124 + 1.5*IT_0127);
    const ccomplex_t IT_0129 = (-0.333333333333333)*IT_0128;
    const ccomplex_t IT_0130 = IT_0121*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0132 = cpow((-2)*s_23 + IT_0017 + IT_0022 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0133 = IT_0130*IT_0131*IT_0132;
    const ccomplex_t IT_0134 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0135 = IT_0040*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0138 = IT_0005*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = m_b*conj(N_d1)*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0141 = IT_0047*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0136 + (-3)*IT_0139 + 3
      *IT_0142);
    const ccomplex_t IT_0144 = 0.166666666666667*IT_0143;
    const ccomplex_t IT_0145 = IT_0052*IT_0144;
    const ccomplex_t IT_0146 = IT_0063*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = m_b*IT_0147;
    const ccomplex_t IT_0149 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0150 = IT_0040*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0153 = IT_0005*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = m_b*conj(N_d2)*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0156 = IT_0047*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0151 + (-3)*IT_0154 + 3
      *IT_0157);
    const ccomplex_t IT_0159 = 0.166666666666667*IT_0158;
    const ccomplex_t IT_0160 = IT_0075*IT_0159;
    const ccomplex_t IT_0161 = IT_0086*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = m_b*IT_0162;
    const ccomplex_t IT_0164 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0165 = IT_0040*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0168 = IT_0005*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = m_b*conj(N_d3)*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0171 = IT_0047*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*(IT_0166 + (-3)*IT_0169 + 3
      *IT_0172);
    const ccomplex_t IT_0174 = 0.166666666666667*IT_0173;
    const ccomplex_t IT_0175 = IT_0098*IT_0174;
    const ccomplex_t IT_0176 = IT_0109*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = m_b*IT_0177;
    const ccomplex_t IT_0179 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0180 = IT_0040*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0183 = IT_0005*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = m_b*conj(N_d4)*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0186 = IT_0047*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0181 + (-3)*IT_0184 + 3
      *IT_0187);
    const ccomplex_t IT_0189 = 0.166666666666667*IT_0188;
    const ccomplex_t IT_0190 = IT_0121*IT_0189;
    const ccomplex_t IT_0191 = IT_0132*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = m_b*IT_0192;
    const ccomplex_t IT_0194 = pow(m_sG, 2);
    const ccomplex_t IT_0195 = cpow((-2)*s_23 + IT_0017 + IT_0022 + -IT_0194 +
       reg_prop, -1);
    const ccomplex_t IT_0196 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0197 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0198 = IT_0196*IT_0197;
    const ccomplex_t IT_0199 = IT_0195*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = m_b*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0203 = IT_0196*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0205 = IT_0195*IT_0203*IT_0204;
    const ccomplex_t IT_0206 = -IT_0064 + -IT_0087 + -IT_0110 + -IT_0133 +
       IT_0148 + IT_0163 + IT_0178 + IT_0193 + (-0.166666666666667)*IT_0201 +
       0.166666666666667*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_10);
    const ccomplex_t IT_0208 = IT_0197*IT_0207;
    const ccomplex_t IT_0209 = IT_0195*IT_0204*IT_0208;
    const ccomplex_t IT_0210 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0211 = IT_0040*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = m_b*conj(N_d1)*e_em*IT_0000*conj(U_sb_00);
    const ccomplex_t IT_0214 = IT_0047*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*(IT_0212 + 1.5*IT_0215);
    const ccomplex_t IT_0217 = (-0.333333333333333)*IT_0216;
    const ccomplex_t IT_0218 = IT_0144*IT_0217;
    const ccomplex_t IT_0219 = IT_0062*IT_0063*IT_0218;
    const ccomplex_t IT_0220 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0221 = IT_0040*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = m_b*conj(N_d2)*e_em*IT_0000*conj(U_sb_00);
    const ccomplex_t IT_0224 = IT_0047*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*(IT_0222 + 1.5*IT_0225);
    const ccomplex_t IT_0227 = (-0.333333333333333)*IT_0226;
    const ccomplex_t IT_0228 = IT_0159*IT_0227;
    const ccomplex_t IT_0229 = IT_0085*IT_0086*IT_0228;
    const ccomplex_t IT_0230 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0231 = IT_0040*IT_0230;
    const ccomplex_t IT_0232 = 1.4142135623731*IT_0231;
    const ccomplex_t IT_0233 = m_b*conj(N_d3)*e_em*IT_0000*conj(U_sb_00);
    const ccomplex_t IT_0234 = IT_0047*IT_0233;
    const ccomplex_t IT_0235 = 1.4142135623731*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*(IT_0232 + 1.5*IT_0235);
    const ccomplex_t IT_0237 = (-0.333333333333333)*IT_0236;
    const ccomplex_t IT_0238 = IT_0174*IT_0237;
    const ccomplex_t IT_0239 = IT_0108*IT_0109*IT_0238;
    const ccomplex_t IT_0240 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0241 = IT_0040*IT_0240;
    const ccomplex_t IT_0242 = 1.4142135623731*IT_0241;
    const ccomplex_t IT_0243 = m_b*conj(N_d4)*e_em*IT_0000*conj(U_sb_00);
    const ccomplex_t IT_0244 = IT_0047*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*(IT_0242 + 1.5*IT_0245);
    const ccomplex_t IT_0247 = (-0.333333333333333)*IT_0246;
    const ccomplex_t IT_0248 = IT_0189*IT_0247;
    const ccomplex_t IT_0249 = IT_0131*IT_0132*IT_0248;
    const ccomplex_t IT_0250 = IT_0060*IT_0217;
    const ccomplex_t IT_0251 = IT_0063*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = m_b*IT_0252;
    const ccomplex_t IT_0254 = IT_0083*IT_0227;
    const ccomplex_t IT_0255 = IT_0086*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = m_b*IT_0256;
    const ccomplex_t IT_0258 = IT_0106*IT_0237;
    const ccomplex_t IT_0259 = IT_0109*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*IT_0259;
    const ccomplex_t IT_0261 = m_b*IT_0260;
    const ccomplex_t IT_0262 = IT_0129*IT_0247;
    const ccomplex_t IT_0263 = IT_0132*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*IT_0263;
    const ccomplex_t IT_0265 = m_b*IT_0264;
    const ccomplex_t IT_0266 = IT_0202*IT_0207;
    const ccomplex_t IT_0267 = IT_0195*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*IT_0267;
    const ccomplex_t IT_0269 = m_b*IT_0268;
    const ccomplex_t IT_0270 = 0.166666666666667*IT_0209 + -IT_0219 + -IT_0229
       + -IT_0239 + -IT_0249 + IT_0253 + IT_0257 + IT_0261 + IT_0265 + (
      -0.166666666666667)*IT_0269;
    const ccomplex_t IT_0271 = IT_0206 + IT_0270;
    const ccomplex_t IT_0272 = -IT_0017;
    const ccomplex_t IT_0273 = s_34 + IT_0272;
    const ccomplex_t IT_0274 = 6*IT_0273;
    const ccomplex_t IT_0275 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0276 = IT_0004*IT_0040;
    const ccomplex_t IT_0277 = e_em*IT_0276;
    const ccomplex_t IT_0278 = IT_0005*IT_0015;
    const ccomplex_t IT_0279 = e_em*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*(IT_0277 + 3*IT_0279);
    const ccomplex_t IT_0281 = (-0.166666666666667)*IT_0280;
    const ccomplex_t IT_0282 = IT_0275*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0005
      *IT_0015 + 1./3*IT_0004*IT_0040)*U_sb_00*conj(U_sb_00) + IT_0004*IT_0040
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0285 = (-0.666666666666667)*IT_0284;
    const ccomplex_t IT_0286 = IT_0283*IT_0285;
    const ccomplex_t IT_0287 = (-0.5)*IT_0200 + -IT_0286;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0040;
    const ccomplex_t IT_0289 = 0.333333333333333*IT_0288;
    const ccomplex_t IT_0290 = IT_0275*IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*IT_0290;
    const ccomplex_t IT_0292 = IT_0285*IT_0291;
    const ccomplex_t IT_0293 = (-0.5)*IT_0268 + -IT_0292;
    const ccomplex_t IT_0294 = IT_0287 + IT_0293;
    const ccomplex_t IT_0295 = -s_14;
    const ccomplex_t IT_0296 = s_13 + IT_0295;
    const ccomplex_t IT_0297 = m_b*IT_0296;
    const ccomplex_t IT_0298 = 18*IT_0297;
    const ccomplex_t IT_0299 = (-0.333333333333333)*IT_0284;
    const ccomplex_t IT_0300 = m_b*IT_0299;
    const ccomplex_t IT_0301 = m_b*IT_0285;
    const ccomplex_t IT_0302 = -IT_0301;
    const ccomplex_t IT_0303 = IT_0300 + IT_0302;
    const ccomplex_t IT_0304 = IT_0283*IT_0303;
    const ccomplex_t IT_0305 = IT_0291*IT_0303;
    const ccomplex_t IT_0306 = (-0.5)*IT_0209 + 0.5*IT_0269 + IT_0304 + 
      -IT_0305;
    const ccomplex_t IT_0307 = (-18)*IT_0017;
    const ccomplex_t IT_0308 = 18*s_34;
    const ccomplex_t IT_0309 = IT_0307 + IT_0308;
    const ccomplex_t IT_0310 = 0.5*IT_0201 + (-0.5)*IT_0205 + -IT_0304 +
       IT_0305;
    const ccomplex_t IT_0311 = m_b*e_em*IT_0000*IT_0002*IT_0005*IT_0011;
    const ccomplex_t IT_0312 = (-0.5)*IT_0311;
    const ccomplex_t IT_0313 = m_b*e_em*mu_h*IT_0000*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0314 = IT_0005*IT_0313;
    const ccomplex_t IT_0315 = m_b*e_em*mu_h*IT_0000*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0316 = IT_0005*IT_0315;
    const ccomplex_t IT_0317 = A_b*IT_0008*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0318 = IT_0011*IT_0317;
    const ccomplex_t IT_0319 = 1.4142135623731*IT_0318;
    const ccomplex_t IT_0320 = A_b*IT_0008*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0321 = IT_0011*IT_0320;
    const ccomplex_t IT_0322 = 1.4142135623731*IT_0321;
    const ccomplex_t IT_0323 = (-0.5)*IT_0314 + 0.5*IT_0316 + (-0.5)*IT_0319 +
       0.5*IT_0322;
    const ccomplex_t IT_0324 = IT_0312*IT_0323;
    const ccomplex_t IT_0325 = cpow((-2)*s_12 + (-2)*IT_0022 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0326 = IT_0324*IT_0325;
    const ccomplex_t IT_0327 = (0 + _Complex_I*1)*IT_0326;
    const ccomplex_t IT_0328 = m_b*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0329 = 0.5*IT_0328;
    const ccomplex_t IT_0330 = IT_0002*IT_0005*IT_0011;
    const ccomplex_t IT_0331 = IT_0313*IT_0330;
    const ccomplex_t IT_0332 = IT_0315*IT_0330;
    const ccomplex_t IT_0333 = IT_0001*IT_0317;
    const ccomplex_t IT_0334 = 1.4142135623731*IT_0333;
    const ccomplex_t IT_0335 = IT_0001*IT_0320;
    const ccomplex_t IT_0336 = 1.4142135623731*IT_0335;
    const ccomplex_t IT_0337 = (-0.5)*IT_0331 + 0.5*IT_0332 + 0.5*IT_0334 + (
      -0.5)*IT_0336;
    const ccomplex_t IT_0338 = IT_0329*IT_0337;
    const ccomplex_t IT_0339 = pow(m_Z, 2);
    const ccomplex_t IT_0340 = cpow((-2)*s_12 + (-2)*IT_0022 + IT_0339 + 
      -reg_prop, -1);
    const ccomplex_t IT_0341 = IT_0338*IT_0340;
    const ccomplex_t IT_0342 = (0 + _Complex_I*1)*IT_0341;
    const ccomplex_t IT_0343 = -IT_0327 + -IT_0342;
    const ccomplex_t IT_0344 = 6*s_34;
    const ccomplex_t IT_0345 = -IT_0252 + -IT_0256 + -IT_0260 + -IT_0264 +
       0.166666666666667*IT_0268;
    const ccomplex_t IT_0346 = m_b*s_13;
    const ccomplex_t IT_0347 = 6*IT_0346;
    const ccomplex_t IT_0348 = 18*IT_0346;
    const ccomplex_t IT_0349 = m_b*s_14;
    const ccomplex_t IT_0350 = (-18)*IT_0349;
    const ccomplex_t IT_0351 = (-6)*IT_0017;
    const ccomplex_t IT_0352 = -IT_0147 + -IT_0162 + -IT_0177 + -IT_0192 +
       0.166666666666667*IT_0200;
    const ccomplex_t IT_0353 = (-6)*IT_0349;
    const ccomplex_t IT_0354 = IT_0306*IT_0308 + IT_0038*IT_0309 + IT_0307
      *IT_0310 + (18*IT_0017 + IT_0308)*IT_0343 + IT_0270*IT_0344 + IT_0345
      *IT_0347 + IT_0293*IT_0348 + IT_0287*IT_0350 + IT_0206*IT_0351 + IT_0352
      *IT_0353;
    const ccomplex_t IT_0355 = (-18)*s_34 + IT_0307;
    const ccomplex_t IT_0356 = conj(IT_0310)*(IT_0306*IT_0307 + IT_0038
      *IT_0309 + IT_0308*IT_0310 + IT_0206*IT_0344 + IT_0287*IT_0348 + IT_0293
      *IT_0350 + IT_0270*IT_0351 + IT_0347*IT_0352 + IT_0345*IT_0353 + IT_0343
      *IT_0355);
    const ccomplex_t IT_0357 = s_34 + IT_0017;
    const ccomplex_t IT_0358 = 6*IT_0357;
    const ccomplex_t IT_0359 = s_13 + s_14;
    const ccomplex_t IT_0360 = m_b*IT_0359;
    const ccomplex_t IT_0361 = 6*IT_0360;
    const ccomplex_t IT_0362 = 18*IT_0360;
    const ccomplex_t IT_0363 = (-18)*IT_0360;
    const ccomplex_t IT_0364 = (-6)*IT_0357;
    const ccomplex_t IT_0365 = (-6)*IT_0360;
    const ccomplex_t IT_0366 = IT_0306*(18*IT_0017 + IT_0308) + (36*IT_0017 +
       IT_0039)*IT_0343 + IT_0310*IT_0355 + IT_0270*IT_0358 + IT_0345*IT_0361 +
       IT_0293*IT_0362 + IT_0287*IT_0363 + IT_0206*IT_0364 + IT_0352*IT_0365;
    const ccomplex_t IT_0367 = cpow(s_12 + IT_0022 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0368 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0369 = IT_0367*IT_0368;
    const ccomplex_t IT_0370 = (0 + _Complex_I*1)*IT_0369;
    const ccomplex_t IT_0371 = (-2)*IT_0368;
    const ccomplex_t IT_0372 = IT_0370*IT_0371;
    const ccomplex_t IT_0373 = 0.25*IT_0372;
    const ccomplex_t IT_0374 = conj(IT_0306) + conj(IT_0310);
    const ccomplex_t IT_0375 = 6*IT_0297;
    const ccomplex_t IT_0376 = conj(IT_0206) + conj(IT_0270);
    const ccomplex_t IT_0377 = (-18)*IT_0297;
    const ccomplex_t IT_0378 = IT_0306 + IT_0310;
    const ccomplex_t IT_0379 = (-6)*IT_0297;
    const ccomplex_t IT_0380 = s_34*IT_0022;
    const ccomplex_t IT_0381 = (-36)*IT_0380;
    const ccomplex_t IT_0382 = IT_0017*IT_0022;
    const ccomplex_t IT_0383 = (-36)*IT_0382;
    const ccomplex_t IT_0384 = s_13*s_14;
    const ccomplex_t IT_0385 = 72*IT_0384;
    const ccomplex_t IT_0386 = IT_0381 + IT_0383 + IT_0385;
    const ccomplex_t IT_0387 = -IT_0386;
    const ccomplex_t IT_0388 = IT_0345 + IT_0352;
    const ccomplex_t IT_0389 = conj(IT_0038)*IT_0375;
    const ccomplex_t IT_0390 = (-18)*IT_0380;
    const ccomplex_t IT_0391 = (-18)*IT_0382;
    const ccomplex_t IT_0392 = 36*IT_0384;
    const ccomplex_t IT_0393 = IT_0390 + IT_0391 + IT_0392;
    const ccomplex_t IT_0394 = -IT_0393;
    const ccomplex_t IT_0395 = IT_0373*IT_0394;
    const ccomplex_t IT_0396 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0397 = 0.666666666666667*IT_0396;
    const ccomplex_t IT_0398 = (-0.333333333333333)*IT_0396;
    const ccomplex_t IT_0399 = IT_0397*IT_0398;
    const ccomplex_t IT_0400 = IT_0367*IT_0399;
    const ccomplex_t IT_0401 = (0 + _Complex_I*1)*IT_0400;
    const ccomplex_t IT_0402 = (-0.0833333333333333)*IT_0372 + 0.5*IT_0401;
    const ccomplex_t IT_0403 = (-6)*IT_0380;
    const ccomplex_t IT_0404 = (-6)*IT_0382;
    const ccomplex_t IT_0405 = 12*IT_0384;
    const ccomplex_t IT_0406 = IT_0403 + IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = -IT_0406;
    const ccomplex_t IT_0408 = ((-36)*IT_0038 + 36*conj(IT_0038))*IT_0297 +
       IT_0298*IT_0374 + IT_0375*IT_0376 + IT_0377*IT_0378 + IT_0271*IT_0379 + (
      -2)*IT_0373*((-12)*IT_0380 + (-12)*IT_0382 + 24*IT_0384) + (conj(IT_0287) 
      + conj(IT_0293))*IT_0393 + IT_0294*IT_0394 + IT_0387*IT_0402 + (conj
      (IT_0345) + conj(IT_0352))*IT_0406 + IT_0388*IT_0407;
    const ccomplex_t IT_0409 = IT_0373*IT_0407;
    const ccomplex_t IT_0410 = conj(IT_0038)*IT_0309;
    const ccomplex_t IT_0411 = IT_0038*IT_0274;
    const ccomplex_t IT_0412 = IT_0403 + IT_0405;
    const ccomplex_t IT_0413 = IT_0390 + IT_0392;
    const ccomplex_t IT_0414 = IT_0373*IT_0393;
    const ccomplex_t IT_0415 = IT_0038*IT_0375;
    const ccomplex_t IT_0416 = conj(IT_0345)*(IT_0306*IT_0347 + IT_0270
      *IT_0348 + IT_0206*IT_0350 + IT_0310*IT_0353 + IT_0343*IT_0361 + IT_0352
      *IT_0391 + IT_0287*IT_0404 + IT_0293*IT_0412 + IT_0345*IT_0413 + IT_0414 +
       IT_0415);
    const ccomplex_t IT_0417 = conj(IT_0352)*(IT_0310*IT_0347 + IT_0206
      *IT_0348 + IT_0270*IT_0350 + IT_0306*IT_0353 + IT_0343*IT_0365 + IT_0345
      *IT_0391 + IT_0293*IT_0404 + IT_0287*IT_0412 + IT_0352*IT_0413 + IT_0414 +
       IT_0415);
    const ccomplex_t IT_0418 = IT_0373*IT_0406;
    const ccomplex_t IT_0419 = IT_0038*IT_0298;
    const ccomplex_t IT_0420 = conj(IT_0293)*(IT_0270*IT_0347 + IT_0306
      *IT_0348 + IT_0310*IT_0350 + IT_0206*IT_0353 + IT_0343*IT_0362 + IT_0287
      *IT_0391 + IT_0352*IT_0404 + IT_0345*IT_0412 + IT_0293*IT_0413 + IT_0418 +
       IT_0419);
    const ccomplex_t IT_0421 = conj(IT_0287)*(IT_0206*IT_0347 + IT_0310
      *IT_0348 + IT_0306*IT_0350 + IT_0270*IT_0353 + IT_0343*IT_0363 + IT_0293
      *IT_0391 + IT_0345*IT_0404 + IT_0352*IT_0412 + IT_0287*IT_0413 + IT_0418 +
       IT_0419);
    const ccomplex_t IT_0422 = conj(IT_0038)*(IT_0038*((-36)*IT_0017 + IT_0039
      ) + IT_0271*IT_0274 + IT_0294*IT_0298) + conj(IT_0306)*IT_0354 + IT_0356 +
       conj(IT_0343)*IT_0366 + IT_0373*(((-12)*IT_0038 + 12*conj(IT_0038))
      *IT_0297 + IT_0374*IT_0375 + IT_0298*IT_0376 + IT_0271*IT_0377 + IT_0378
      *IT_0379 + IT_0373*IT_0387) + IT_0388*(IT_0389 + IT_0395) + IT_0402
      *IT_0408 + IT_0294*IT_0409 + IT_0378*IT_0410 + conj(IT_0270)*(IT_0206
      *IT_0307 + IT_0270*IT_0308 + IT_0306*IT_0344 + IT_0293*IT_0347 + IT_0345
      *IT_0348 + IT_0310*IT_0351 + IT_0350*IT_0352 + IT_0287*IT_0353 + IT_0343
      *IT_0358 + IT_0411) + conj(IT_0206)*(IT_0270*IT_0307 + IT_0206*IT_0308 +
       IT_0310*IT_0344 + IT_0287*IT_0347 + IT_0345*IT_0350 + IT_0306*IT_0351 +
       IT_0348*IT_0352 + IT_0293*IT_0353 + IT_0343*IT_0364 + IT_0411) + IT_0416 
      + IT_0417 + IT_0420 + IT_0421;
    return create_ccomplex_return(IT_0422);
}

