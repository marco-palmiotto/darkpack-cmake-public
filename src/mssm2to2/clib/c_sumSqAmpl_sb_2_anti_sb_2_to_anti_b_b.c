#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_anti_sb_2_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_anti_sb_2_to_anti_b_b(
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
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
    const ccomplex_t IT_0018 = (-2)*U_sb_01*(((-0.5)*m_W*IT_0014*(IT_0005 + 1.
      /3*IT_0004*IT_0016) + IT_0000*IT_0002*IT_0003*IT_0005*IT_0017)*conj
      (U_sb_01) + -1./2*IT_0002*IT_0005*IT_0009*IT_0010*conj(U_sb_11)) + U_sb_11
      *(IT_0002*IT_0005*IT_0009*IT_0010*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0004*IT_0014*IT_0016 + (-3)*IT_0000*IT_0002*IT_0003*IT_0005*IT_0017)
      *conj(U_sb_11));
    const ccomplex_t IT_0019 = (0 + _Complex_I*-1)*(e_em*IT_0018 + (
      -1.4142135623731)*A_b*IT_0003*IT_0008*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = pow(m_sb_2, 2);
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
    const ccomplex_t IT_0031 = 2*U_sb_01*((IT_0000*IT_0002*IT_0005*IT_0010
      *IT_0017 + (-0.5)*m_W*(IT_0005 + 1./3*IT_0004*IT_0016)*IT_0030)*conj
      (U_sb_01) + 1./2*IT_0002*IT_0003*IT_0005*IT_0009*conj(U_sb_11)) + U_sb_11*
      (IT_0002*IT_0003*IT_0005*IT_0009*conj(U_sb_01) + 2*(IT_0000*IT_0002
      *IT_0005*IT_0010*IT_0017 + -1./3*m_W*IT_0004*IT_0016*IT_0030)*conj(U_sb_11
      ));
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(e_em*IT_0031 +
       1.4142135623731*A_b*IT_0008*IT_0010*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -IT_0025 + -IT_0037;
    const ccomplex_t IT_0039 = pow(m_sG, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0017 + IT_0022 + -IT_0039 +
       reg_prop, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_11);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = IT_0040*IT_0041*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_01);
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = IT_0040*IT_0041*IT_0049;
    const ccomplex_t IT_0051 = cpow(IT_0015, -1);
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0005
      *IT_0015 + 1./3*IT_0004*IT_0051)*U_sb_01*conj(U_sb_01) + IT_0004*IT_0051
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0053 = (-0.333333333333333)*IT_0052;
    const ccomplex_t IT_0054 = m_b*IT_0053;
    const ccomplex_t IT_0055 = (-0.666666666666667)*IT_0052;
    const ccomplex_t IT_0056 = m_b*IT_0055;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = IT_0054 + IT_0057;
    const ccomplex_t IT_0059 = IT_0004*IT_0051;
    const ccomplex_t IT_0060 = e_em*IT_0059;
    const ccomplex_t IT_0061 = IT_0005*IT_0015;
    const ccomplex_t IT_0062 = e_em*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0060 + 3*IT_0062);
    const ccomplex_t IT_0064 = (-0.166666666666667)*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0058*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0051;
    const ccomplex_t IT_0070 = 0.333333333333333*IT_0069;
    const ccomplex_t IT_0071 = IT_0065*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0058*IT_0072;
    const ccomplex_t IT_0074 = IT_0043*IT_0047;
    const ccomplex_t IT_0075 = IT_0040*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = m_b*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0050 + -IT_0068 + IT_0073 + 0.5
      *IT_0077;
    const ccomplex_t IT_0079 = IT_0042*IT_0048;
    const ccomplex_t IT_0080 = IT_0040*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = m_b*IT_0081;
    const ccomplex_t IT_0083 = IT_0068 + -IT_0073 + 0.5*IT_0082;
    const ccomplex_t IT_0084 = 18*s_34;
    const ccomplex_t IT_0085 = (-18)*IT_0017;
    const ccomplex_t IT_0086 = IT_0084 + IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0088 = (-2)*IT_0087;
    const ccomplex_t IT_0089 = cpow(s_12 + IT_0022 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0088*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0094 = 0.666666666666667*IT_0093;
    const ccomplex_t IT_0095 = (-0.333333333333333)*IT_0093;
    const ccomplex_t IT_0096 = IT_0094*IT_0095;
    const ccomplex_t IT_0097 = IT_0089*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = (-0.0833333333333333)*IT_0092 + 0.5*IT_0098;
    const ccomplex_t IT_0100 = cpow(IT_0099, 2);
    const ccomplex_t IT_0101 = s_34*IT_0022;
    const ccomplex_t IT_0102 = (-36)*IT_0101;
    const ccomplex_t IT_0103 = IT_0017*IT_0022;
    const ccomplex_t IT_0104 = (-36)*IT_0103;
    const ccomplex_t IT_0105 = s_13*s_14;
    const ccomplex_t IT_0106 = 72*IT_0105;
    const ccomplex_t IT_0107 = IT_0102 + IT_0104 + IT_0106;
    const ccomplex_t IT_0108 = 0.25*IT_0092;
    const ccomplex_t IT_0109 = (-6)*IT_0101;
    const ccomplex_t IT_0110 = 12*IT_0105;
    const ccomplex_t IT_0111 = IT_0109 + IT_0110;
    const ccomplex_t IT_0112 = (-0.5)*IT_0076;
    const ccomplex_t IT_0113 = IT_0055*IT_0067;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = conj(IT_0112) + conj(IT_0114);
    const ccomplex_t IT_0116 = cpow((-2)*s_23 + IT_0017 + IT_0022 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0117 = N_B1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0118 = IT_0051*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = N_W1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0121 = IT_0005*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = IT_0002*IT_0005;
    const ccomplex_t IT_0124 = m_b*N_d1*e_em*IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0119 + (-3)*IT_0122 + 3
      *IT_0126);
    const ccomplex_t IT_0128 = 0.166666666666667*IT_0127;
    const ccomplex_t IT_0129 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0130 = IT_0051*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0133 = IT_0005*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = m_b*conj(N_d1)*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0136 = IT_0123*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0131 + (-3)*IT_0134 + 3
      *IT_0137);
    const ccomplex_t IT_0139 = 0.166666666666667*IT_0138;
    const ccomplex_t IT_0140 = IT_0128*IT_0139;
    const ccomplex_t IT_0141 = IT_0116*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = cpow((-2)*s_23 + IT_0017 + IT_0022 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0144 = N_B2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0145 = IT_0051*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = N_W2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0148 = IT_0005*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = m_b*N_d2*e_em*IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0151 = IT_0123*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0146 + (-3)*IT_0149 + 3
      *IT_0152);
    const ccomplex_t IT_0154 = 0.166666666666667*IT_0153;
    const ccomplex_t IT_0155 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0156 = IT_0051*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0159 = IT_0005*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = m_b*conj(N_d2)*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0162 = IT_0123*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0157 + (-3)*IT_0160 + 3
      *IT_0163);
    const ccomplex_t IT_0165 = 0.166666666666667*IT_0164;
    const ccomplex_t IT_0166 = IT_0154*IT_0165;
    const ccomplex_t IT_0167 = IT_0143*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = cpow((-2)*s_23 + IT_0017 + IT_0022 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0170 = N_B4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0171 = IT_0051*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = N_W4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0174 = IT_0005*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = m_b*N_d4*e_em*IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0177 = IT_0123*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*(IT_0172 + (-3)*IT_0175 + 3
      *IT_0178);
    const ccomplex_t IT_0180 = 0.166666666666667*IT_0179;
    const ccomplex_t IT_0181 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0182 = IT_0051*IT_0181;
    const ccomplex_t IT_0183 = 1.4142135623731*IT_0182;
    const ccomplex_t IT_0184 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0185 = IT_0005*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = m_b*conj(N_d4)*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0188 = IT_0123*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*(IT_0183 + (-3)*IT_0186 + 3
      *IT_0189);
    const ccomplex_t IT_0191 = 0.166666666666667*IT_0190;
    const ccomplex_t IT_0192 = IT_0180*IT_0191;
    const ccomplex_t IT_0193 = IT_0169*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = cpow((-2)*s_23 + IT_0017 + IT_0022 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0196 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0197 = IT_0051*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0200 = IT_0005*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = m_b*N_d3*e_em*IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0203 = IT_0123*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*(IT_0198 + (-3)*IT_0201 + 3
      *IT_0204);
    const ccomplex_t IT_0206 = 0.166666666666667*IT_0205;
    const ccomplex_t IT_0207 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0208 = IT_0051*IT_0207;
    const ccomplex_t IT_0209 = 1.4142135623731*IT_0208;
    const ccomplex_t IT_0210 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0211 = IT_0005*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = m_b*conj(N_d3)*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0214 = IT_0123*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*(IT_0209 + (-3)*IT_0212 + 3
      *IT_0215);
    const ccomplex_t IT_0217 = 0.166666666666667*IT_0216;
    const ccomplex_t IT_0218 = IT_0206*IT_0217;
    const ccomplex_t IT_0219 = IT_0195*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = -IT_0142 + -IT_0168 + -IT_0194 + -IT_0220;
    const ccomplex_t IT_0222 = 0.166666666666667*IT_0076;
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = IT_0112 + IT_0114;
    const ccomplex_t IT_0225 = conj(IT_0221) + conj(IT_0222);
    const ccomplex_t IT_0226 = (-0.5)*IT_0081;
    const ccomplex_t IT_0227 = IT_0055*IT_0072;
    const ccomplex_t IT_0228 = -IT_0227;
    const ccomplex_t IT_0229 = conj(IT_0226) + conj(IT_0228);
    const ccomplex_t IT_0230 = conj(N_B4)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0231 = IT_0051*IT_0230;
    const ccomplex_t IT_0232 = 1.4142135623731*IT_0231;
    const ccomplex_t IT_0233 = m_b*conj(N_d4)*e_em*IT_0000*conj(U_sb_01);
    const ccomplex_t IT_0234 = IT_0123*IT_0233;
    const ccomplex_t IT_0235 = 1.4142135623731*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*(IT_0232 + 1.5*IT_0235);
    const ccomplex_t IT_0237 = (-0.333333333333333)*IT_0236;
    const ccomplex_t IT_0238 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0239 = IT_0051*IT_0238;
    const ccomplex_t IT_0240 = 1.4142135623731*IT_0239;
    const ccomplex_t IT_0241 = m_b*N_d4*e_em*IT_0000*U_sb_01;
    const ccomplex_t IT_0242 = IT_0123*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*(IT_0240 + 1.5*IT_0243);
    const ccomplex_t IT_0245 = (-0.333333333333333)*IT_0244;
    const ccomplex_t IT_0246 = IT_0237*IT_0245;
    const ccomplex_t IT_0247 = IT_0169*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = conj(N_B2)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0250 = IT_0051*IT_0249;
    const ccomplex_t IT_0251 = 1.4142135623731*IT_0250;
    const ccomplex_t IT_0252 = m_b*conj(N_d2)*e_em*IT_0000*conj(U_sb_01);
    const ccomplex_t IT_0253 = IT_0123*IT_0252;
    const ccomplex_t IT_0254 = 1.4142135623731*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*(IT_0251 + 1.5*IT_0254);
    const ccomplex_t IT_0256 = (-0.333333333333333)*IT_0255;
    const ccomplex_t IT_0257 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0258 = IT_0051*IT_0257;
    const ccomplex_t IT_0259 = 1.4142135623731*IT_0258;
    const ccomplex_t IT_0260 = m_b*N_d2*e_em*IT_0000*U_sb_01;
    const ccomplex_t IT_0261 = IT_0123*IT_0260;
    const ccomplex_t IT_0262 = 1.4142135623731*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*(IT_0259 + 1.5*IT_0262);
    const ccomplex_t IT_0264 = (-0.333333333333333)*IT_0263;
    const ccomplex_t IT_0265 = IT_0256*IT_0264;
    const ccomplex_t IT_0266 = IT_0143*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0269 = IT_0051*IT_0268;
    const ccomplex_t IT_0270 = 1.4142135623731*IT_0269;
    const ccomplex_t IT_0271 = m_b*conj(N_d3)*e_em*IT_0000*conj(U_sb_01);
    const ccomplex_t IT_0272 = IT_0123*IT_0271;
    const ccomplex_t IT_0273 = 1.4142135623731*IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*(IT_0270 + 1.5*IT_0273);
    const ccomplex_t IT_0275 = (-0.333333333333333)*IT_0274;
    const ccomplex_t IT_0276 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0277 = IT_0051*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = m_b*N_d3*e_em*IT_0000*U_sb_01;
    const ccomplex_t IT_0280 = IT_0123*IT_0279;
    const ccomplex_t IT_0281 = 1.4142135623731*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*(IT_0278 + 1.5*IT_0281);
    const ccomplex_t IT_0283 = (-0.333333333333333)*IT_0282;
    const ccomplex_t IT_0284 = IT_0275*IT_0283;
    const ccomplex_t IT_0285 = IT_0195*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = -IT_0248 + -IT_0267 + -IT_0286;
    const ccomplex_t IT_0288 = conj(N_B1)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0289 = IT_0051*IT_0288;
    const ccomplex_t IT_0290 = 1.4142135623731*IT_0289;
    const ccomplex_t IT_0291 = m_b*conj(N_d1)*e_em*IT_0000*conj(U_sb_01);
    const ccomplex_t IT_0292 = IT_0123*IT_0291;
    const ccomplex_t IT_0293 = 1.4142135623731*IT_0292;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*(IT_0290 + 1.5*IT_0293);
    const ccomplex_t IT_0295 = (-0.333333333333333)*IT_0294;
    const ccomplex_t IT_0296 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0297 = IT_0051*IT_0296;
    const ccomplex_t IT_0298 = 1.4142135623731*IT_0297;
    const ccomplex_t IT_0299 = m_b*N_d1*e_em*IT_0000*U_sb_01;
    const ccomplex_t IT_0300 = IT_0123*IT_0299;
    const ccomplex_t IT_0301 = 1.4142135623731*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*(IT_0298 + 1.5*IT_0301);
    const ccomplex_t IT_0303 = (-0.333333333333333)*IT_0302;
    const ccomplex_t IT_0304 = IT_0295*IT_0303;
    const ccomplex_t IT_0305 = IT_0116*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*IT_0305;
    const ccomplex_t IT_0307 = 0.166666666666667*IT_0081 + -IT_0306;
    const ccomplex_t IT_0308 = IT_0287 + IT_0307;
    const ccomplex_t IT_0309 = IT_0226 + IT_0228;
    const ccomplex_t IT_0310 = conj(IT_0287) + conj(IT_0307);
    const ccomplex_t IT_0311 = (-18)*IT_0101;
    const ccomplex_t IT_0312 = 36*IT_0105;
    const ccomplex_t IT_0313 = IT_0311 + IT_0312;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0315 = IT_0191*IT_0237;
    const ccomplex_t IT_0316 = IT_0169*IT_0314*IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0318 = IT_0139*IT_0295;
    const ccomplex_t IT_0319 = IT_0116*IT_0317*IT_0318;
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0321 = IT_0165*IT_0256;
    const ccomplex_t IT_0322 = IT_0143*IT_0320*IT_0321;
    const ccomplex_t IT_0323 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0324 = IT_0217*IT_0275;
    const ccomplex_t IT_0325 = IT_0195*IT_0323*IT_0324;
    const ccomplex_t IT_0326 = m_b*IT_0267;
    const ccomplex_t IT_0327 = m_b*IT_0306;
    const ccomplex_t IT_0328 = m_b*IT_0286;
    const ccomplex_t IT_0329 = m_b*IT_0248;
    const ccomplex_t IT_0330 = 0.166666666666667*IT_0045 + (-0.166666666666667
      )*IT_0082 + -IT_0316 + -IT_0319 + -IT_0322 + -IT_0325 + IT_0326 + IT_0327 
      + IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = IT_0206*IT_0283;
    const ccomplex_t IT_0332 = IT_0195*IT_0323*IT_0331;
    const ccomplex_t IT_0333 = IT_0180*IT_0245;
    const ccomplex_t IT_0334 = IT_0169*IT_0314*IT_0333;
    const ccomplex_t IT_0335 = IT_0128*IT_0303;
    const ccomplex_t IT_0336 = IT_0116*IT_0317*IT_0335;
    const ccomplex_t IT_0337 = IT_0154*IT_0264;
    const ccomplex_t IT_0338 = IT_0143*IT_0320*IT_0337;
    const ccomplex_t IT_0339 = m_b*IT_0142;
    const ccomplex_t IT_0340 = m_b*IT_0168;
    const ccomplex_t IT_0341 = m_b*IT_0220;
    const ccomplex_t IT_0342 = m_b*IT_0194;
    const ccomplex_t IT_0343 = 0.166666666666667*IT_0050 + (-0.166666666666667
      )*IT_0077 + -IT_0332 + -IT_0334 + -IT_0336 + -IT_0338 + IT_0339 + IT_0340 
      + IT_0341 + IT_0342;
    const ccomplex_t IT_0344 = (-18)*IT_0343;
    const ccomplex_t IT_0345 = (-6)*IT_0078;
    const ccomplex_t IT_0346 = IT_0344 + IT_0345;
    const ccomplex_t IT_0347 = (-18)*conj(IT_0343);
    const ccomplex_t IT_0348 = (-6)*conj(IT_0078);
    const ccomplex_t IT_0349 = IT_0347 + IT_0348;
    const ccomplex_t IT_0350 = conj(IT_0046) + conj(IT_0083);
    const ccomplex_t IT_0351 = (-18)*IT_0078;
    const ccomplex_t IT_0352 = (-6)*IT_0343;
    const ccomplex_t IT_0353 = IT_0351 + IT_0352;
    const ccomplex_t IT_0354 = IT_0046 + IT_0083;
    const ccomplex_t IT_0355 = (-18)*conj(IT_0078);
    const ccomplex_t IT_0356 = (-6)*conj(IT_0343);
    const ccomplex_t IT_0357 = IT_0355 + IT_0356;
    const ccomplex_t IT_0358 = 6*s_34;
    const ccomplex_t IT_0359 = -IT_0017;
    const ccomplex_t IT_0360 = 6*IT_0038;
    const ccomplex_t IT_0361 = 6*conj(IT_0038);
    const ccomplex_t IT_0362 = m_b*s_13;
    const ccomplex_t IT_0363 = 6*IT_0078;
    const ccomplex_t IT_0364 = 18*IT_0343;
    const ccomplex_t IT_0365 = 6*conj(IT_0078);
    const ccomplex_t IT_0366 = 18*conj(IT_0343);
    const ccomplex_t IT_0367 = 6*IT_0343;
    const ccomplex_t IT_0368 = 18*IT_0078;
    const ccomplex_t IT_0369 = 6*conj(IT_0343);
    const ccomplex_t IT_0370 = 18*conj(IT_0078);
    const ccomplex_t IT_0371 = 6*IT_0046;
    const ccomplex_t IT_0372 = 6*IT_0083;
    const ccomplex_t IT_0373 = 18*IT_0330;
    const ccomplex_t IT_0374 = 6*conj(IT_0046);
    const ccomplex_t IT_0375 = 6*conj(IT_0083);
    const ccomplex_t IT_0376 = 18*conj(IT_0330);
    const ccomplex_t IT_0377 = 6*IT_0330;
    const ccomplex_t IT_0378 = 18*IT_0046;
    const ccomplex_t IT_0379 = 18*IT_0083;
    const ccomplex_t IT_0380 = 6*conj(IT_0330);
    const ccomplex_t IT_0381 = 18*conj(IT_0046);
    const ccomplex_t IT_0382 = 18*conj(IT_0083);
    const ccomplex_t IT_0383 = IT_0225*(IT_0363 + IT_0364) + IT_0223*(IT_0365 
      + IT_0366) + IT_0115*(IT_0367 + IT_0368) + IT_0224*(IT_0369 + IT_0370) +
       IT_0310*(IT_0371 + IT_0372 + IT_0373) + IT_0308*(IT_0374 + IT_0375 +
       IT_0376) + IT_0229*(IT_0377 + IT_0378 + IT_0379) + IT_0309*(IT_0380 +
       IT_0381 + IT_0382);
    const ccomplex_t IT_0384 = m_b*e_em*IT_0000*IT_0002*IT_0005*IT_0011;
    const ccomplex_t IT_0385 = (-0.5)*IT_0384;
    const ccomplex_t IT_0386 = m_b*e_em*mu_h*IT_0000*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0387 = IT_0005*IT_0386;
    const ccomplex_t IT_0388 = m_b*e_em*mu_h*IT_0000*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0389 = IT_0005*IT_0388;
    const ccomplex_t IT_0390 = A_b*IT_0008*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0391 = IT_0011*IT_0390;
    const ccomplex_t IT_0392 = 1.4142135623731*IT_0391;
    const ccomplex_t IT_0393 = A_b*IT_0008*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0394 = IT_0011*IT_0393;
    const ccomplex_t IT_0395 = 1.4142135623731*IT_0394;
    const ccomplex_t IT_0396 = (-0.5)*IT_0387 + 0.5*IT_0389 + (-0.5)*IT_0392 +
       0.5*IT_0395;
    const ccomplex_t IT_0397 = IT_0385*IT_0396;
    const ccomplex_t IT_0398 = cpow((-2)*s_12 + (-2)*IT_0022 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0399 = IT_0397*IT_0398;
    const ccomplex_t IT_0400 = (0 + _Complex_I*1)*IT_0399;
    const ccomplex_t IT_0401 = pow(m_Z, 2);
    const ccomplex_t IT_0402 = cpow((-2)*s_12 + (-2)*IT_0022 + IT_0401 + 
      -reg_prop, -1);
    const ccomplex_t IT_0403 = m_b*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0404 = 0.5*IT_0403;
    const ccomplex_t IT_0405 = IT_0002*IT_0005*IT_0011;
    const ccomplex_t IT_0406 = IT_0386*IT_0405;
    const ccomplex_t IT_0407 = IT_0388*IT_0405;
    const ccomplex_t IT_0408 = IT_0001*IT_0390;
    const ccomplex_t IT_0409 = 1.4142135623731*IT_0408;
    const ccomplex_t IT_0410 = IT_0001*IT_0393;
    const ccomplex_t IT_0411 = 1.4142135623731*IT_0410;
    const ccomplex_t IT_0412 = (-0.5)*IT_0406 + 0.5*IT_0407 + 0.5*IT_0409 + (
      -0.5)*IT_0411;
    const ccomplex_t IT_0413 = IT_0404*IT_0412;
    const ccomplex_t IT_0414 = IT_0402*IT_0413;
    const ccomplex_t IT_0415 = (0 + _Complex_I*1)*IT_0414;
    const ccomplex_t IT_0416 = -IT_0400 + -IT_0415;
    const ccomplex_t IT_0417 = (-18)*s_34 + IT_0085;
    const ccomplex_t IT_0418 = 18*IT_0017;
    const ccomplex_t IT_0419 = IT_0084 + IT_0418;
    const ccomplex_t IT_0420 = s_34 + IT_0017;
    const ccomplex_t IT_0421 = -s_14;
    const ccomplex_t IT_0422 = s_13 + IT_0421;
    const ccomplex_t IT_0423 = m_b*IT_0422;
    const ccomplex_t IT_0424 = IT_0225 + conj(IT_0287) + conj(IT_0307);
    const ccomplex_t IT_0425 = IT_0223 + IT_0287 + IT_0307;
    const ccomplex_t IT_0426 = IT_0115 + conj(IT_0226) + conj(IT_0228);
    const ccomplex_t IT_0427 = 18*IT_0038;
    const ccomplex_t IT_0428 = IT_0224 + IT_0226 + IT_0228;
    const ccomplex_t IT_0429 = 18*conj(IT_0038);
    const ccomplex_t IT_0430 = (-18)*IT_0046;
    const ccomplex_t IT_0431 = (-18)*IT_0083;
    const ccomplex_t IT_0432 = (-6)*IT_0330;
    const ccomplex_t IT_0433 = (-18)*IT_0330;
    const ccomplex_t IT_0434 = (-6)*IT_0046;
    const ccomplex_t IT_0435 = (-6)*IT_0083;
    const ccomplex_t IT_0436 = IT_0360*IT_0424 + IT_0361*IT_0425 + IT_0426
      *IT_0427 + IT_0428*IT_0429 + IT_0099*((-36)*IT_0038 + 36*conj(IT_0038) +
       IT_0351 + IT_0352 + IT_0369 + IT_0370 + IT_0380 + IT_0381 + IT_0382 +
       IT_0430 + IT_0431 + IT_0432) + IT_0108*((-12)*IT_0038 + 12*conj(IT_0038) 
      + IT_0344 + IT_0345 + IT_0365 + IT_0366 + IT_0374 + IT_0375 + IT_0376 +
       IT_0433 + IT_0434 + IT_0435);
    const ccomplex_t IT_0437 = m_b*s_14;
    const ccomplex_t IT_0438 = (-18)*conj(IT_0046);
    const ccomplex_t IT_0439 = (-18)*conj(IT_0083);
    const ccomplex_t IT_0440 = (-6)*conj(IT_0330);
    const ccomplex_t IT_0441 = (-18)*conj(IT_0330);
    const ccomplex_t IT_0442 = (-6)*conj(IT_0046);
    const ccomplex_t IT_0443 = (-6)*conj(IT_0083);
    const ccomplex_t IT_0444 = IT_0310*IT_0346 + IT_0308*IT_0349 + IT_0229
      *IT_0353 + IT_0309*IT_0357 + IT_0115*(IT_0430 + IT_0431 + IT_0432) +
       IT_0225*(IT_0433 + IT_0434 + IT_0435) + IT_0224*(IT_0438 + IT_0439 +
       IT_0440) + IT_0223*(IT_0441 + IT_0442 + IT_0443);
    const ccomplex_t IT_0445 = s_13 + s_14;
    const ccomplex_t IT_0446 = m_b*IT_0445;
    const ccomplex_t IT_0447 = 6*IT_0416;
    const ccomplex_t IT_0448 = 6*conj(IT_0416);
    const ccomplex_t IT_0449 = 18*IT_0416;
    const ccomplex_t IT_0450 = 18*conj(IT_0416);
    const ccomplex_t IT_0451 = (-18)*IT_0416;
    const ccomplex_t IT_0452 = (-18)*conj(IT_0416);
    const ccomplex_t IT_0453 = (-6)*IT_0416;
    const ccomplex_t IT_0454 = (-6)*conj(IT_0416);
    const ccomplex_t IT_0455 = (-18)*IT_0221;
    const ccomplex_t IT_0456 = (-18)*IT_0222;
    const ccomplex_t IT_0457 = (-6)*IT_0112;
    const ccomplex_t IT_0458 = (-6)*IT_0114;
    const ccomplex_t IT_0459 = (-18)*conj(IT_0221);
    const ccomplex_t IT_0460 = (-18)*conj(IT_0222);
    const ccomplex_t IT_0461 = (-6)*conj(IT_0112);
    const ccomplex_t IT_0462 = (-6)*conj(IT_0114);
    const ccomplex_t IT_0463 = (-18)*IT_0112;
    const ccomplex_t IT_0464 = (-18)*IT_0114;
    const ccomplex_t IT_0465 = (-6)*IT_0221;
    const ccomplex_t IT_0466 = (-6)*IT_0222;
    const ccomplex_t IT_0467 = (-18)*conj(IT_0112);
    const ccomplex_t IT_0468 = (-18)*conj(IT_0114);
    const ccomplex_t IT_0469 = (-6)*conj(IT_0221);
    const ccomplex_t IT_0470 = (-6)*conj(IT_0222);
    const ccomplex_t IT_0471 = IT_0103*(IT_0310*(IT_0455 + IT_0456 + IT_0457 +
       IT_0458) + IT_0308*(IT_0459 + IT_0460 + IT_0461 + IT_0462) + IT_0229*
      (IT_0463 + IT_0464 + IT_0465 + IT_0466) + IT_0309*(IT_0467 + IT_0468 +
       IT_0469 + IT_0470));
    const ccomplex_t IT_0472 = (-18)*IT_0103;
    const ccomplex_t IT_0473 = IT_0311 + IT_0312 + IT_0472;
    const ccomplex_t IT_0474 = -IT_0108;
    const ccomplex_t IT_0475 = -IT_0099;
    const ccomplex_t IT_0476 = (-6)*IT_0103;
    const ccomplex_t IT_0477 = IT_0109 + IT_0110 + IT_0476;
    const ccomplex_t IT_0478 = 36*s_34;
    const ccomplex_t IT_0479 = (conj(IT_0038)*(IT_0046 + IT_0078 + IT_0083) +
       IT_0038*(conj(IT_0046) + conj(IT_0078) + conj(IT_0083)))*IT_0086 + 
      -IT_0100*IT_0107 + (-2)*IT_0108*(IT_0099*((-12)*IT_0101 + (-12)*IT_0103 +
       24*IT_0105) + 0.5*IT_0107*IT_0108) + IT_0111*(IT_0115*IT_0223 + IT_0224
      *IT_0225 + IT_0229*IT_0308 + IT_0309*IT_0310) + (IT_0115*IT_0224 + IT_0223
      *IT_0225 + IT_0229*IT_0309 + IT_0308*IT_0310)*IT_0313 + IT_0084*(IT_0078
      *conj(IT_0078) + IT_0330*conj(IT_0330) + IT_0343*conj(IT_0343) + (IT_0046 
      + IT_0083)*IT_0350) + IT_0017*(conj(IT_0330)*IT_0346 + IT_0330*IT_0349 +
       IT_0350*IT_0353 + IT_0354*IT_0357) + (conj(IT_0078)*IT_0343 + IT_0078
      *conj(IT_0343) + IT_0330*IT_0350 + conj(IT_0330)*IT_0354)*IT_0358 + (s_34 
      + IT_0359)*((conj(IT_0330) + conj(IT_0343))*IT_0360 + (IT_0330 + IT_0343)
      *IT_0361) + IT_0362*IT_0383 + conj(IT_0416)*(IT_0078*IT_0417 + IT_0354
      *IT_0419 + (IT_0352 + IT_0377)*IT_0420) + IT_0423*IT_0436 + IT_0437
      *IT_0444 + IT_0446*(IT_0310*IT_0447 + IT_0308*IT_0448 + IT_0229*IT_0449 +
       IT_0309*IT_0450 + IT_0115*IT_0451 + IT_0224*IT_0452 + IT_0225*IT_0453 +
       IT_0223*IT_0454) + IT_0471 + IT_0473*(IT_0108*IT_0424 + IT_0099*IT_0426 +
       IT_0425*IT_0474 + IT_0428*IT_0475) + (IT_0099*IT_0424 + IT_0108*IT_0426 +
       IT_0428*IT_0474 + IT_0425*IT_0475)*IT_0477 + IT_0038*conj(IT_0038)*((-36)
      *IT_0017 + IT_0478) + IT_0416*(conj(IT_0078)*IT_0417 + IT_0350*IT_0419 + 
      (IT_0356 + IT_0380)*IT_0420 + conj(IT_0416)*(36*IT_0017 + IT_0478));
    return create_ccomplex_return(IT_0479);
}

