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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_01*conj(U_sb_01) + IT_0001*IT_0002
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0005 = (-0.333333333333333)*IT_0004;
    const ccomplex_t IT_0006 = m_b*IT_0005;
    const ccomplex_t IT_0007 = (-0.666666666666667)*IT_0004;
    const ccomplex_t IT_0008 = m_b*IT_0007;
    const ccomplex_t IT_0009 = -IT_0008;
    const ccomplex_t IT_0010 = IT_0006 + IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0012 = 0.333333333333333*IT_0011;
    const ccomplex_t IT_0013 = pow(m_sb_2, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = pow(m_b, 2);
    const ccomplex_t IT_0019 = pow(m_sG, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0013 + IT_0018 + -IT_0019 +
       reg_prop, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_01);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = m_b*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0028 = IT_0021*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0030 = IT_0020*IT_0028*IT_0029;
    const ccomplex_t IT_0031 = IT_0001*IT_0002;
    const ccomplex_t IT_0032 = e_em*IT_0031;
    const ccomplex_t IT_0033 = IT_0000*IT_0003;
    const ccomplex_t IT_0034 = e_em*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0032 + 3*IT_0034);
    const ccomplex_t IT_0036 = (-0.166666666666667)*IT_0035;
    const ccomplex_t IT_0037 = IT_0014*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0010*IT_0038;
    const ccomplex_t IT_0040 = IT_0017 + 0.5*IT_0026 + (-0.5)*IT_0030 + 
      -IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_11);
    const ccomplex_t IT_0042 = IT_0027*IT_0041;
    const ccomplex_t IT_0043 = IT_0020*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = m_b*IT_0044;
    const ccomplex_t IT_0046 = IT_0022*IT_0041;
    const ccomplex_t IT_0047 = IT_0020*IT_0029*IT_0046;
    const ccomplex_t IT_0048 = -IT_0017 + IT_0039 + 0.5*IT_0045 + (-0.5)
      *IT_0047;
    const ccomplex_t IT_0049 = pow(m_W, -1);
    const ccomplex_t IT_0050 = cos(beta);
    const ccomplex_t IT_0051 = cpow(IT_0050, -1);
    const ccomplex_t IT_0052 = sin(alpha);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_b*e_em*IT_0003*IT_0049
      *IT_0051*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = cpow(V_tb, 2);
    const ccomplex_t IT_0056 = m_b*mu_h*IT_0049;
    const ccomplex_t IT_0057 = cos(alpha);
    const ccomplex_t IT_0058 = IT_0050*IT_0052;
    const ccomplex_t IT_0059 = sin(beta);
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = IT_0058 + IT_0060;
    const ccomplex_t IT_0062 = cpow(IT_0000, -2);
    const ccomplex_t IT_0063 = 2*U_sb_01*((IT_0003*IT_0018*IT_0049*IT_0051
      *IT_0052 + (-0.5)*m_W*IT_0061*(IT_0003 + 1./3*IT_0002*IT_0062))*conj
      (U_sb_01) + 1./2*IT_0003*IT_0051*IT_0056*IT_0057*conj(U_sb_11)) + U_sb_11*
      (IT_0003*IT_0051*IT_0056*IT_0057*conj(U_sb_01) + 2*(IT_0003*IT_0018
      *IT_0049*IT_0051*IT_0052 + -1./3*m_W*IT_0002*IT_0061*IT_0062)*conj(U_sb_11
      ));
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(e_em*IT_0063 +
       1.4142135623731*A_b*IT_0052*IT_0055*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0054*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_b*e_em*IT_0003*IT_0049
      *IT_0051*IT_0057;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = -IT_0052*IT_0059;
    const ccomplex_t IT_0073 = IT_0050*IT_0057;
    const ccomplex_t IT_0074 = IT_0072 + IT_0073;
    const ccomplex_t IT_0075 = (-2)*U_sb_01*((IT_0003*IT_0018*IT_0049*IT_0051
      *IT_0057 + (-0.5)*m_W*(IT_0003 + 1./3*IT_0002*IT_0062)*IT_0074)*conj
      (U_sb_01) + -1./2*IT_0003*IT_0051*IT_0052*IT_0056*conj(U_sb_11)) + U_sb_11
      *(IT_0003*IT_0051*IT_0052*IT_0056*conj(U_sb_01) + (-2)*(IT_0003*IT_0018
      *IT_0049*IT_0051*IT_0057 + -1./3*m_W*IT_0002*IT_0062*IT_0074)*conj(U_sb_11
      ));
    const ccomplex_t IT_0076 = (0 + _Complex_I*-1)*(e_em*IT_0075 + (
      -1.4142135623731)*A_b*IT_0055*IT_0057*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = IT_0071*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + (-2)*IT_0013 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = -IT_0069 + -IT_0081;
    const ccomplex_t IT_0083 = 18*s_34;
    const ccomplex_t IT_0084 = (-18)*IT_0018;
    const ccomplex_t IT_0085 = IT_0083 + IT_0084;
    const ccomplex_t IT_0086 = IT_0007*IT_0038;
    const ccomplex_t IT_0087 = (-0.5)*IT_0025 + -IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_23 + IT_0013 + IT_0018 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0089 = N_B1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0090 = IT_0001*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = N_W1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0093 = IT_0003*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = IT_0003*IT_0051;
    const ccomplex_t IT_0096 = m_b*N_d1*e_em*IT_0049*conj(U_sb_11);
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0091 + (-3)*IT_0094 + 3
      *IT_0098);
    const ccomplex_t IT_0100 = 0.166666666666667*IT_0099;
    const ccomplex_t IT_0101 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0102 = IT_0001*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0105 = IT_0003*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = m_b*conj(N_d1)*e_em*IT_0049*U_sb_11;
    const ccomplex_t IT_0108 = IT_0095*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + (-3)*IT_0106 + 3
      *IT_0109);
    const ccomplex_t IT_0111 = 0.166666666666667*IT_0110;
    const ccomplex_t IT_0112 = IT_0100*IT_0111;
    const ccomplex_t IT_0113 = IT_0088*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_23 + IT_0013 + IT_0018 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0116 = N_B2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0117 = IT_0001*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = N_W2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0120 = IT_0003*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = m_b*N_d2*e_em*IT_0049*conj(U_sb_11);
    const ccomplex_t IT_0123 = IT_0095*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0118 + (-3)*IT_0121 + 3
      *IT_0124);
    const ccomplex_t IT_0126 = 0.166666666666667*IT_0125;
    const ccomplex_t IT_0127 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0128 = IT_0001*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0131 = IT_0003*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = m_b*conj(N_d2)*e_em*IT_0049*U_sb_11;
    const ccomplex_t IT_0134 = IT_0095*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0129 + (-3)*IT_0132 + 3
      *IT_0135);
    const ccomplex_t IT_0137 = 0.166666666666667*IT_0136;
    const ccomplex_t IT_0138 = IT_0126*IT_0137;
    const ccomplex_t IT_0139 = IT_0115*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = cpow((-2)*s_23 + IT_0013 + IT_0018 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0142 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0143 = IT_0001*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0146 = IT_0003*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = m_b*N_d3*e_em*IT_0049*conj(U_sb_11);
    const ccomplex_t IT_0149 = IT_0095*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0144 + (-3)*IT_0147 + 3
      *IT_0150);
    const ccomplex_t IT_0152 = 0.166666666666667*IT_0151;
    const ccomplex_t IT_0153 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0154 = IT_0001*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0157 = IT_0003*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = m_b*conj(N_d3)*e_em*IT_0049*U_sb_11;
    const ccomplex_t IT_0160 = IT_0095*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*(IT_0155 + (-3)*IT_0158 + 3
      *IT_0161);
    const ccomplex_t IT_0163 = 0.166666666666667*IT_0162;
    const ccomplex_t IT_0164 = IT_0152*IT_0163;
    const ccomplex_t IT_0165 = IT_0141*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = N_B4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0168 = IT_0001*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = N_W4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0171 = IT_0003*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = m_b*N_d4*e_em*IT_0049*conj(U_sb_11);
    const ccomplex_t IT_0174 = IT_0095*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0169 + (-3)*IT_0172 + 3
      *IT_0175);
    const ccomplex_t IT_0177 = 0.166666666666667*IT_0176;
    const ccomplex_t IT_0178 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0179 = IT_0001*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0182 = IT_0003*IT_0181;
    const ccomplex_t IT_0183 = 1.4142135623731*IT_0182;
    const ccomplex_t IT_0184 = m_b*conj(N_d4)*e_em*IT_0049*U_sb_11;
    const ccomplex_t IT_0185 = IT_0095*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0180 + (-3)*IT_0183 + 3
      *IT_0186);
    const ccomplex_t IT_0188 = 0.166666666666667*IT_0187;
    const ccomplex_t IT_0189 = IT_0177*IT_0188;
    const ccomplex_t IT_0190 = cpow((-2)*s_23 + IT_0013 + IT_0018 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0191 = IT_0189*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = 0.166666666666667*IT_0025 + -IT_0114 + -IT_0140
       + -IT_0166 + -IT_0192;
    const ccomplex_t IT_0194 = IT_0007*IT_0016;
    const ccomplex_t IT_0195 = (-0.5)*IT_0044 + -IT_0194;
    const ccomplex_t IT_0196 = conj(N_B2)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0197 = IT_0001*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = m_b*conj(N_d2)*e_em*IT_0049*conj(U_sb_01);
    const ccomplex_t IT_0200 = IT_0095*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*(IT_0198 + 1.5*IT_0201);
    const ccomplex_t IT_0203 = (-0.333333333333333)*IT_0202;
    const ccomplex_t IT_0204 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0205 = IT_0001*IT_0204;
    const ccomplex_t IT_0206 = 1.4142135623731*IT_0205;
    const ccomplex_t IT_0207 = m_b*N_d2*e_em*IT_0049*U_sb_01;
    const ccomplex_t IT_0208 = IT_0095*IT_0207;
    const ccomplex_t IT_0209 = 1.4142135623731*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*(IT_0206 + 1.5*IT_0209);
    const ccomplex_t IT_0211 = (-0.333333333333333)*IT_0210;
    const ccomplex_t IT_0212 = IT_0203*IT_0211;
    const ccomplex_t IT_0213 = IT_0115*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0216 = IT_0001*IT_0215;
    const ccomplex_t IT_0217 = 1.4142135623731*IT_0216;
    const ccomplex_t IT_0218 = m_b*conj(N_d3)*e_em*IT_0049*conj(U_sb_01);
    const ccomplex_t IT_0219 = IT_0095*IT_0218;
    const ccomplex_t IT_0220 = 1.4142135623731*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*(IT_0217 + 1.5*IT_0220);
    const ccomplex_t IT_0222 = (-0.333333333333333)*IT_0221;
    const ccomplex_t IT_0223 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0224 = IT_0001*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = m_b*N_d3*e_em*IT_0049*U_sb_01;
    const ccomplex_t IT_0227 = IT_0095*IT_0226;
    const ccomplex_t IT_0228 = 1.4142135623731*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*(IT_0225 + 1.5*IT_0228);
    const ccomplex_t IT_0230 = (-0.333333333333333)*IT_0229;
    const ccomplex_t IT_0231 = IT_0222*IT_0230;
    const ccomplex_t IT_0232 = IT_0141*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = conj(N_B1)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0235 = IT_0001*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = m_b*conj(N_d1)*e_em*IT_0049*conj(U_sb_01);
    const ccomplex_t IT_0238 = IT_0095*IT_0237;
    const ccomplex_t IT_0239 = 1.4142135623731*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*(IT_0236 + 1.5*IT_0239);
    const ccomplex_t IT_0241 = (-0.333333333333333)*IT_0240;
    const ccomplex_t IT_0242 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0243 = IT_0001*IT_0242;
    const ccomplex_t IT_0244 = 1.4142135623731*IT_0243;
    const ccomplex_t IT_0245 = m_b*N_d1*e_em*IT_0049*U_sb_01;
    const ccomplex_t IT_0246 = IT_0095*IT_0245;
    const ccomplex_t IT_0247 = 1.4142135623731*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*(IT_0244 + 1.5*IT_0247);
    const ccomplex_t IT_0249 = (-0.333333333333333)*IT_0248;
    const ccomplex_t IT_0250 = IT_0241*IT_0249;
    const ccomplex_t IT_0251 = IT_0088*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = conj(N_B4)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0254 = IT_0001*IT_0253;
    const ccomplex_t IT_0255 = 1.4142135623731*IT_0254;
    const ccomplex_t IT_0256 = m_b*conj(N_d4)*e_em*IT_0049*conj(U_sb_01);
    const ccomplex_t IT_0257 = IT_0095*IT_0256;
    const ccomplex_t IT_0258 = 1.4142135623731*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*(IT_0255 + 1.5*IT_0258);
    const ccomplex_t IT_0260 = (-0.333333333333333)*IT_0259;
    const ccomplex_t IT_0261 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0262 = IT_0001*IT_0261;
    const ccomplex_t IT_0263 = 1.4142135623731*IT_0262;
    const ccomplex_t IT_0264 = m_b*N_d4*e_em*IT_0049*U_sb_01;
    const ccomplex_t IT_0265 = IT_0095*IT_0264;
    const ccomplex_t IT_0266 = 1.4142135623731*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*(IT_0263 + 1.5*IT_0266);
    const ccomplex_t IT_0268 = (-0.333333333333333)*IT_0267;
    const ccomplex_t IT_0269 = IT_0260*IT_0268;
    const ccomplex_t IT_0270 = IT_0190*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = 0.166666666666667*IT_0044 + -IT_0214 + -IT_0233
       + -IT_0252 + -IT_0271;
    const ccomplex_t IT_0273 = s_34*IT_0013;
    const ccomplex_t IT_0274 = (-6)*IT_0273;
    const ccomplex_t IT_0275 = s_13*s_14;
    const ccomplex_t IT_0276 = 12*IT_0275;
    const ccomplex_t IT_0277 = IT_0274 + IT_0276;
    const ccomplex_t IT_0278 = (-18)*IT_0273;
    const ccomplex_t IT_0279 = 36*IT_0275;
    const ccomplex_t IT_0280 = IT_0278 + IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0282 = (-2)*IT_0281;
    const ccomplex_t IT_0283 = cpow(s_12 + IT_0013 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0284 = IT_0281*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = IT_0282*IT_0285;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0288 = 0.666666666666667*IT_0287;
    const ccomplex_t IT_0289 = (-0.333333333333333)*IT_0287;
    const ccomplex_t IT_0290 = IT_0288*IT_0289;
    const ccomplex_t IT_0291 = IT_0283*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = (-0.0833333333333333)*IT_0286 + 0.5*IT_0292;
    const ccomplex_t IT_0294 = cpow(IT_0293, 2);
    const ccomplex_t IT_0295 = (-36)*IT_0273;
    const ccomplex_t IT_0296 = IT_0013*IT_0018;
    const ccomplex_t IT_0297 = (-36)*IT_0296;
    const ccomplex_t IT_0298 = 72*IT_0275;
    const ccomplex_t IT_0299 = IT_0295 + IT_0297 + IT_0298;
    const ccomplex_t IT_0300 = 0.25*IT_0286;
    const ccomplex_t IT_0301 = (-18)*IT_0040;
    const ccomplex_t IT_0302 = IT_0177*IT_0268;
    const ccomplex_t IT_0303 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0304 = IT_0190*IT_0302*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0306 = IT_0100*IT_0249;
    const ccomplex_t IT_0307 = IT_0088*IT_0305*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0309 = IT_0126*IT_0211;
    const ccomplex_t IT_0310 = IT_0115*IT_0308*IT_0309;
    const ccomplex_t IT_0311 = IT_0152*IT_0230;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0313 = IT_0141*IT_0311*IT_0312;
    const ccomplex_t IT_0314 = m_b*IT_0114;
    const ccomplex_t IT_0315 = m_b*IT_0140;
    const ccomplex_t IT_0316 = m_b*IT_0166;
    const ccomplex_t IT_0317 = m_b*IT_0192;
    const ccomplex_t IT_0318 = (-0.166666666666667)*IT_0026 +
       0.166666666666667*IT_0030 + -IT_0304 + -IT_0307 + -IT_0310 + -IT_0313 +
       IT_0314 + IT_0315 + IT_0316 + IT_0317;
    const ccomplex_t IT_0319 = (-6)*IT_0318;
    const ccomplex_t IT_0320 = IT_0301 + IT_0319;
    const ccomplex_t IT_0321 = (-18)*conj(IT_0040);
    const ccomplex_t IT_0322 = (-6)*conj(IT_0318);
    const ccomplex_t IT_0323 = IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = IT_0188*IT_0260;
    const ccomplex_t IT_0325 = IT_0190*IT_0303*IT_0324;
    const ccomplex_t IT_0326 = IT_0111*IT_0241;
    const ccomplex_t IT_0327 = IT_0088*IT_0305*IT_0326;
    const ccomplex_t IT_0328 = IT_0137*IT_0203;
    const ccomplex_t IT_0329 = IT_0115*IT_0308*IT_0328;
    const ccomplex_t IT_0330 = IT_0163*IT_0222;
    const ccomplex_t IT_0331 = IT_0141*IT_0312*IT_0330;
    const ccomplex_t IT_0332 = m_b*IT_0252;
    const ccomplex_t IT_0333 = m_b*IT_0214;
    const ccomplex_t IT_0334 = m_b*IT_0233;
    const ccomplex_t IT_0335 = m_b*IT_0271;
    const ccomplex_t IT_0336 = (-0.166666666666667)*IT_0045 +
       0.166666666666667*IT_0047 + -IT_0325 + -IT_0327 + -IT_0329 + -IT_0331 +
       IT_0332 + IT_0333 + IT_0334 + IT_0335;
    const ccomplex_t IT_0337 = (-18)*IT_0318;
    const ccomplex_t IT_0338 = (-6)*IT_0040;
    const ccomplex_t IT_0339 = IT_0337 + IT_0338;
    const ccomplex_t IT_0340 = (-18)*conj(IT_0318);
    const ccomplex_t IT_0341 = (-6)*conj(IT_0040);
    const ccomplex_t IT_0342 = IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = 6*s_34;
    const ccomplex_t IT_0344 = -IT_0018;
    const ccomplex_t IT_0345 = 6*IT_0082;
    const ccomplex_t IT_0346 = 6*conj(IT_0082);
    const ccomplex_t IT_0347 = m_b*e_em*IT_0003*IT_0049*IT_0051*IT_0059;
    const ccomplex_t IT_0348 = (-0.5)*IT_0347;
    const ccomplex_t IT_0349 = m_b*e_em*mu_h*IT_0049*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0350 = IT_0003*IT_0349;
    const ccomplex_t IT_0351 = m_b*e_em*mu_h*IT_0049*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0352 = IT_0003*IT_0351;
    const ccomplex_t IT_0353 = A_b*IT_0055*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0354 = IT_0059*IT_0353;
    const ccomplex_t IT_0355 = 1.4142135623731*IT_0354;
    const ccomplex_t IT_0356 = A_b*IT_0055*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0357 = IT_0059*IT_0356;
    const ccomplex_t IT_0358 = 1.4142135623731*IT_0357;
    const ccomplex_t IT_0359 = (-0.5)*IT_0350 + 0.5*IT_0352 + (-0.5)*IT_0355 +
       0.5*IT_0358;
    const ccomplex_t IT_0360 = IT_0348*IT_0359;
    const ccomplex_t IT_0361 = cpow((-2)*s_12 + (-2)*IT_0013 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0362 = IT_0360*IT_0361;
    const ccomplex_t IT_0363 = (0 + _Complex_I*1)*IT_0362;
    const ccomplex_t IT_0364 = m_b*e_em*IT_0003*IT_0049;
    const ccomplex_t IT_0365 = 0.5*IT_0364;
    const ccomplex_t IT_0366 = IT_0003*IT_0051*IT_0059;
    const ccomplex_t IT_0367 = IT_0349*IT_0366;
    const ccomplex_t IT_0368 = IT_0351*IT_0366;
    const ccomplex_t IT_0369 = IT_0050*IT_0353;
    const ccomplex_t IT_0370 = 1.4142135623731*IT_0369;
    const ccomplex_t IT_0371 = IT_0050*IT_0356;
    const ccomplex_t IT_0372 = 1.4142135623731*IT_0371;
    const ccomplex_t IT_0373 = (-0.5)*IT_0367 + 0.5*IT_0368 + 0.5*IT_0370 + (
      -0.5)*IT_0372;
    const ccomplex_t IT_0374 = IT_0365*IT_0373;
    const ccomplex_t IT_0375 = pow(m_Z, 2);
    const ccomplex_t IT_0376 = cpow((-2)*s_12 + (-2)*IT_0013 + IT_0375 + 
      -reg_prop, -1);
    const ccomplex_t IT_0377 = IT_0374*IT_0376;
    const ccomplex_t IT_0378 = (0 + _Complex_I*1)*IT_0377;
    const ccomplex_t IT_0379 = -IT_0363 + -IT_0378;
    const ccomplex_t IT_0380 = 18*IT_0018;
    const ccomplex_t IT_0381 = IT_0083 + IT_0380;
    const ccomplex_t IT_0382 = (-18)*s_34 + IT_0084;
    const ccomplex_t IT_0383 = s_13 + s_14;
    const ccomplex_t IT_0384 = m_b*IT_0383;
    const ccomplex_t IT_0385 = (-6)*IT_0193;
    const ccomplex_t IT_0386 = (-18)*IT_0087 + IT_0385;
    const ccomplex_t IT_0387 = s_34 + IT_0018;
    const ccomplex_t IT_0388 = 6*IT_0336;
    const ccomplex_t IT_0389 = m_b*s_14;
    const ccomplex_t IT_0390 = (-18)*IT_0336;
    const ccomplex_t IT_0391 = (-6)*IT_0048;
    const ccomplex_t IT_0392 = (-18)*IT_0048;
    const ccomplex_t IT_0393 = (-6)*IT_0336;
    const ccomplex_t IT_0394 = (-0.166666666666667)*conj(IT_0195)*IT_0320 + (
      -0.166666666666667)*IT_0195*IT_0323 + (IT_0087 + 3*IT_0193)*conj(IT_0336) 
      + (-0.166666666666667)*conj(IT_0272)*IT_0339 + (-0.166666666666667)
      *IT_0272*IT_0342 + (-0.166666666666667)*conj(IT_0048)*IT_0386 + (
      -0.166666666666667)*conj(IT_0193)*(IT_0390 + IT_0391) + (
      -0.166666666666667)*conj(IT_0087)*(IT_0392 + IT_0393);
    const ccomplex_t IT_0395 = -s_14;
    const ccomplex_t IT_0396 = s_13 + IT_0395;
    const ccomplex_t IT_0397 = m_b*IT_0396;
    const ccomplex_t IT_0398 = conj(IT_0193) + conj(IT_0272);
    const ccomplex_t IT_0399 = IT_0193 + IT_0272;
    const ccomplex_t IT_0400 = conj(IT_0087) + conj(IT_0195);
    const ccomplex_t IT_0401 = 18*IT_0082;
    const ccomplex_t IT_0402 = IT_0087 + IT_0195;
    const ccomplex_t IT_0403 = 18*conj(IT_0082);
    const ccomplex_t IT_0404 = 6*conj(IT_0040);
    const ccomplex_t IT_0405 = 6*conj(IT_0048);
    const ccomplex_t IT_0406 = 18*conj(IT_0318);
    const ccomplex_t IT_0407 = 18*conj(IT_0336);
    const ccomplex_t IT_0408 = 6*conj(IT_0318);
    const ccomplex_t IT_0409 = 6*conj(IT_0336);
    const ccomplex_t IT_0410 = 18*conj(IT_0040);
    const ccomplex_t IT_0411 = 18*conj(IT_0048);
    const ccomplex_t IT_0412 = IT_0345*IT_0398 + IT_0346*IT_0399 + IT_0400
      *IT_0401 + IT_0402*IT_0403 + IT_0300*((-12)*IT_0082 + 12*conj(IT_0082) +
       IT_0337 + IT_0338 + IT_0390 + IT_0391 + IT_0404 + IT_0405 + IT_0406 +
       IT_0407) + IT_0293*((-36)*IT_0082 + 36*conj(IT_0082) + IT_0301 + IT_0319 
      + IT_0392 + IT_0393 + IT_0408 + IT_0409 + IT_0410 + IT_0411);
    const ccomplex_t IT_0413 = (-18)*IT_0296;
    const ccomplex_t IT_0414 = IT_0278 + IT_0279 + IT_0413;
    const ccomplex_t IT_0415 = -IT_0300;
    const ccomplex_t IT_0416 = -IT_0293;
    const ccomplex_t IT_0417 = (-6)*IT_0296;
    const ccomplex_t IT_0418 = IT_0274 + IT_0276 + IT_0417;
    const ccomplex_t IT_0419 = 36*s_34;
    const ccomplex_t IT_0420 = m_b*s_13;
    const ccomplex_t IT_0421 = 6*conj(IT_0272);
    const ccomplex_t IT_0422 = 18*conj(IT_0195);
    const ccomplex_t IT_0423 = 0.166666666666667*IT_0040*(18*conj(IT_0087) + 6
      *conj(IT_0193)) + (conj(IT_0087) + 3*conj(IT_0193))*IT_0318 + 3*conj
      (IT_0272)*IT_0336 + 0.166666666666667*conj(IT_0195)*IT_0388 +
       0.166666666666667*IT_0193*(IT_0404 + IT_0406) + 0.166666666666667*IT_0272
      *(IT_0405 + IT_0407) + 0.166666666666667*IT_0087*(IT_0408 + IT_0410) +
       0.166666666666667*IT_0195*(IT_0409 + IT_0411) + 0.166666666666667*IT_0048
      *(IT_0421 + IT_0422);
    const ccomplex_t IT_0424 = (-6)*conj(IT_0193);
    const ccomplex_t IT_0425 = ((conj(IT_0040) + conj(IT_0048))*IT_0082 + 
      (IT_0040 + IT_0048)*conj(IT_0082))*IT_0085 + (conj(IT_0087)*IT_0193 +
       IT_0087*conj(IT_0193) + conj(IT_0195)*IT_0272 + IT_0195*conj(IT_0272))
      *IT_0277 + (IT_0087*conj(IT_0087) + IT_0193*conj(IT_0193) + IT_0195*conj
      (IT_0195) + IT_0272*conj(IT_0272))*IT_0280 + -IT_0294*IT_0299 + (-2)
      *IT_0300*(IT_0293*((-12)*IT_0273 + 24*IT_0275 + (-12)*IT_0296) + 0.5
      *IT_0299*IT_0300) + IT_0083*(IT_0040*conj(IT_0040) + IT_0048*conj(IT_0048)
       + IT_0318*conj(IT_0318) + IT_0336*conj(IT_0336)) + IT_0018*(conj(IT_0048)
      *IT_0320 + IT_0048*IT_0323 + conj(IT_0336)*IT_0339 + IT_0336*IT_0342) + 
      (conj(IT_0040)*IT_0318 + IT_0040*conj(IT_0318) + conj(IT_0048)*IT_0336 +
       IT_0048*conj(IT_0336))*IT_0343 + (s_34 + IT_0344)*((conj(IT_0318) + conj
      (IT_0336))*IT_0345 + (IT_0318 + IT_0336)*IT_0346) + conj(IT_0379)*(IT_0048
      *IT_0381 + IT_0040*IT_0382 + IT_0384*(18*IT_0195 + 6*IT_0272 + IT_0386) +
       IT_0387*(IT_0319 + IT_0388)) + (-6)*IT_0389*IT_0394 + IT_0397*IT_0412 +
       IT_0414*(IT_0300*IT_0398 + IT_0293*IT_0400 + IT_0399*IT_0415 + IT_0402
      *IT_0416) + (IT_0293*IT_0398 + IT_0300*IT_0400 + IT_0402*IT_0415 + IT_0399
      *IT_0416)*IT_0418 + IT_0082*conj(IT_0082)*((-36)*IT_0018 + IT_0419) + 6
      *IT_0420*IT_0423 + IT_0379*(conj(IT_0048)*IT_0381 + conj(IT_0040)*IT_0382 
      + IT_0387*(IT_0322 + IT_0409) + conj(IT_0379)*(36*IT_0018 + IT_0419) +
       IT_0384*((-18)*conj(IT_0087) + IT_0421 + IT_0422 + IT_0424)) + IT_0296*
      (conj(IT_0087)*((-18)*IT_0195 + (-6)*IT_0272) + IT_0087*((-18)*conj
      (IT_0195) + (-6)*conj(IT_0272)) + (conj(IT_0195) + 3*conj(IT_0272))
      *IT_0385 + (IT_0195 + 3*IT_0272)*IT_0424);
    return create_ccomplex_return(IT_0425);
}

