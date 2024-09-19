#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_sb_1_to_b_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_sb_1_to_b_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
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
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0013 = IT_0003*IT_0012;
    const ccomplex_t IT_0014 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = IT_0005 + IT_0011 + IT_0013 + IT_0015;
    const ccomplex_t IT_0017 = IT_0002*IT_0006;
    const ccomplex_t IT_0018 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0000*IT_0008;
    const ccomplex_t IT_0021 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0022 + -IT_0024 + -IT_0026;
    const ccomplex_t IT_0028 = IT_0016 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = m_N_1*IT_0030;
    const ccomplex_t IT_0032 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0035 = IT_0009*IT_0034;
    const ccomplex_t IT_0036 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0039 = IT_0009*IT_0038;
    const ccomplex_t IT_0040 = IT_0033 + IT_0035 + IT_0037 + IT_0039;
    const ccomplex_t IT_0041 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0042 = IT_0017*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0044 = IT_0020*IT_0043;
    const ccomplex_t IT_0045 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0046 = IT_0017*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0048 = IT_0020*IT_0047;
    const ccomplex_t IT_0049 = -IT_0042 + -IT_0044 + -IT_0046 + -IT_0048;
    const ccomplex_t IT_0050 = IT_0040 + IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = m_N_2*IT_0052;
    const ccomplex_t IT_0054 = IT_0031 + IT_0053;
    const ccomplex_t IT_0055 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0056 = IT_0002*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = cos(beta);
    const ccomplex_t IT_0059 = cpow(IT_0058, -1);
    const ccomplex_t IT_0060 = IT_0008*IT_0059;
    const ccomplex_t IT_0061 = pow(m_W, -1);
    const ccomplex_t IT_0062 = m_b*N_d1*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0057 + 1.5*IT_0064);
    const ccomplex_t IT_0066 = (-0.333333333333333)*IT_0065;
    const ccomplex_t IT_0067 = pow(m_b, 2);
    const ccomplex_t IT_0068 = pow(m_sb_1, 2);
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0067 + IT_0068 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0066*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = IT_0054*IT_0071;
    const ccomplex_t IT_0073 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0074 = IT_0003*IT_0073;
    const ccomplex_t IT_0075 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0076 = IT_0009*IT_0075;
    const ccomplex_t IT_0077 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0078 = IT_0003*IT_0077;
    const ccomplex_t IT_0079 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0080 = IT_0009*IT_0079;
    const ccomplex_t IT_0081 = IT_0074 + IT_0076 + IT_0078 + IT_0080;
    const ccomplex_t IT_0082 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0083 = IT_0017*IT_0082;
    const ccomplex_t IT_0084 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0085 = IT_0020*IT_0084;
    const ccomplex_t IT_0086 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0087 = IT_0017*IT_0086;
    const ccomplex_t IT_0088 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0089 = IT_0020*IT_0088;
    const ccomplex_t IT_0090 = -IT_0083 + -IT_0085 + -IT_0087 + -IT_0089;
    const ccomplex_t IT_0091 = IT_0081 + IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = m_N_3*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0096 = IT_0003*IT_0095;
    const ccomplex_t IT_0097 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0098 = IT_0009*IT_0097;
    const ccomplex_t IT_0099 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0100 = IT_0003*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0102 = IT_0009*IT_0101;
    const ccomplex_t IT_0103 = IT_0096 + IT_0098 + IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0105 = IT_0017*IT_0104;
    const ccomplex_t IT_0106 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0107 = IT_0020*IT_0106;
    const ccomplex_t IT_0108 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0109 = IT_0017*IT_0108;
    const ccomplex_t IT_0110 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0111 = IT_0020*IT_0110;
    const ccomplex_t IT_0112 = -IT_0105 + -IT_0107 + -IT_0109 + -IT_0111;
    const ccomplex_t IT_0113 = IT_0103 + IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = m_N_2*IT_0115;
    const ccomplex_t IT_0117 = IT_0094 + IT_0116;
    const ccomplex_t IT_0118 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0119 = IT_0002*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = m_b*N_d3*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0122 = IT_0060*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + 1.5*IT_0123);
    const ccomplex_t IT_0125 = (-0.333333333333333)*IT_0124;
    const ccomplex_t IT_0126 = cpow((-2)*s_23 + IT_0067 + IT_0068 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0117*IT_0128;
    const ccomplex_t IT_0130 = pow(m_N_2, 2);
    const ccomplex_t IT_0131 = cpow((-2)*s_13 + IT_0067 + IT_0130 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0132 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0133 = IT_0002*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = m_b*N_d2*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0136 = IT_0060*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0134 + 1.5*IT_0137);
    const ccomplex_t IT_0139 = (-0.333333333333333)*IT_0138;
    const ccomplex_t IT_0140 = cpow(V_tb, 2);
    const ccomplex_t IT_0141 = m_b*mu_h*IT_0061;
    const ccomplex_t IT_0142 = sin(beta);
    const ccomplex_t IT_0143 = -IT_0006*IT_0142;
    const ccomplex_t IT_0144 = IT_0000*IT_0058;
    const ccomplex_t IT_0145 = IT_0143 + IT_0144;
    const ccomplex_t IT_0146 = cpow(IT_0001, -2);
    const ccomplex_t IT_0147 = (-2)*U_sb_00*((IT_0000*IT_0008*IT_0059*IT_0061
      *IT_0067 + (-0.5)*m_W*IT_0145*(IT_0008 + 1./3*IT_0007*IT_0146))*conj
      (U_sb_00) + -1./2*IT_0006*IT_0008*IT_0059*IT_0141*conj(U_sb_10)) + U_sb_10
      *(IT_0006*IT_0008*IT_0059*IT_0141*conj(U_sb_00) + (-2)*(IT_0000*IT_0008
      *IT_0059*IT_0061*IT_0067 + -1./3*m_W*IT_0007*IT_0145*IT_0146)*conj(U_sb_10
      ));
    const ccomplex_t IT_0148 = (0 + _Complex_I*-1)*(e_em*IT_0147 + (
      -1.4142135623731)*A_b*IT_0000*IT_0140*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = IT_0139*IT_0149;
    const ccomplex_t IT_0151 = IT_0131*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = cpow((-2)*s_13 + IT_0067 + IT_0130 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0154 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0155 = IT_0002*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = m_b*N_d2*e_em*IT_0061*U_sb_01;
    const ccomplex_t IT_0158 = IT_0060*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*(IT_0156 + 1.5*IT_0159);
    const ccomplex_t IT_0161 = (-0.333333333333333)*IT_0160;
    const ccomplex_t IT_0162 = (-2)*U_sb_00*((IT_0000*IT_0008*IT_0059*IT_0061
      *IT_0067 + (-0.5)*m_W*IT_0145*(IT_0008 + 1./3*IT_0007*IT_0146))*conj
      (U_sb_01) + -1./2*IT_0006*IT_0008*IT_0059*IT_0141*conj(U_sb_11)) + U_sb_10
      *(IT_0006*IT_0008*IT_0059*IT_0141*conj(U_sb_01) + (-2)*(IT_0000*IT_0008
      *IT_0059*IT_0061*IT_0067 + -1./3*m_W*IT_0007*IT_0145*IT_0146)*conj(U_sb_11
      ));
    const ccomplex_t IT_0163 = (0 + _Complex_I*-1)*(e_em*IT_0162 + (
      -1.4142135623731)*A_b*IT_0000*IT_0140*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = IT_0161*IT_0164;
    const ccomplex_t IT_0166 = IT_0153*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = cpow((-2)*s_12 + IT_0067 + -IT_0068 + -IT_0130 
      + -reg_prop, -1);
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0008
      *IT_0059*IT_0061;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = IT_0139*IT_0170;
    const ccomplex_t IT_0172 = IT_0168*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = m_b*IT_0173;
    const ccomplex_t IT_0175 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0176 = IT_0003*IT_0175;
    const ccomplex_t IT_0177 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0178 = IT_0017*IT_0177;
    const ccomplex_t IT_0179 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0180 = IT_0020*IT_0179;
    const ccomplex_t IT_0181 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0182 = IT_0009*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0176 + -IT_0178 + 
      -IT_0180 + IT_0182);
    const ccomplex_t IT_0184 = m_N_2*IT_0183;
    const ccomplex_t IT_0185 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0186 = IT_0003*IT_0185;
    const ccomplex_t IT_0187 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0188 = IT_0017*IT_0187;
    const ccomplex_t IT_0189 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0190 = IT_0020*IT_0189;
    const ccomplex_t IT_0191 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0192 = IT_0009*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*(IT_0186 + -IT_0188 + 
      -IT_0190 + IT_0192);
    const ccomplex_t IT_0194 = m_N_2*IT_0193;
    const ccomplex_t IT_0195 = IT_0184 + IT_0194;
    const ccomplex_t IT_0196 = cpow((-2)*s_23 + IT_0067 + IT_0068 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0197 = IT_0139*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = IT_0195*IT_0198;
    const ccomplex_t IT_0200 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0201 = IT_0002*IT_0200;
    const ccomplex_t IT_0202 = 1.4142135623731*IT_0201;
    const ccomplex_t IT_0203 = m_b*N_d4*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0204 = IT_0060*IT_0203;
    const ccomplex_t IT_0205 = 1.4142135623731*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*(IT_0202 + 1.5*IT_0205);
    const ccomplex_t IT_0207 = (-0.333333333333333)*IT_0206;
    const ccomplex_t IT_0208 = cpow((-2)*s_23 + IT_0067 + IT_0068 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0209 = IT_0207*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0212 = IT_0003*IT_0211;
    const ccomplex_t IT_0213 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0214 = IT_0009*IT_0213;
    const ccomplex_t IT_0215 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0216 = IT_0003*IT_0215;
    const ccomplex_t IT_0217 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0218 = IT_0009*IT_0217;
    const ccomplex_t IT_0219 = IT_0212 + IT_0214 + IT_0216 + IT_0218;
    const ccomplex_t IT_0220 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0221 = IT_0017*IT_0220;
    const ccomplex_t IT_0222 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0223 = IT_0020*IT_0222;
    const ccomplex_t IT_0224 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0225 = IT_0017*IT_0224;
    const ccomplex_t IT_0226 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0227 = IT_0020*IT_0226;
    const ccomplex_t IT_0228 = -IT_0221 + -IT_0223 + -IT_0225 + -IT_0227;
    const ccomplex_t IT_0229 = IT_0219 + IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = 0.5*IT_0230;
    const ccomplex_t IT_0232 = m_N_4*IT_0231;
    const ccomplex_t IT_0233 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0234 = IT_0003*IT_0233;
    const ccomplex_t IT_0235 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0236 = IT_0009*IT_0235;
    const ccomplex_t IT_0237 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0238 = IT_0003*IT_0237;
    const ccomplex_t IT_0239 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0240 = IT_0009*IT_0239;
    const ccomplex_t IT_0241 = IT_0234 + IT_0236 + IT_0238 + IT_0240;
    const ccomplex_t IT_0242 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0243 = IT_0017*IT_0242;
    const ccomplex_t IT_0244 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0245 = IT_0020*IT_0244;
    const ccomplex_t IT_0246 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0247 = IT_0017*IT_0246;
    const ccomplex_t IT_0248 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0249 = IT_0020*IT_0248;
    const ccomplex_t IT_0250 = -IT_0243 + -IT_0245 + -IT_0247 + -IT_0249;
    const ccomplex_t IT_0251 = IT_0241 + IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = 0.5*IT_0252;
    const ccomplex_t IT_0254 = m_N_2*IT_0253;
    const ccomplex_t IT_0255 = IT_0232 + IT_0254;
    const ccomplex_t IT_0256 = IT_0210*IT_0255;
    const ccomplex_t IT_0257 = IT_0072 + IT_0129 + IT_0152 + IT_0167 + (-2)
      *IT_0174 + IT_0199 + IT_0256;
    const ccomplex_t IT_0258 = m_b*conj(N_d2)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0259 = IT_0060*IT_0258;
    const ccomplex_t IT_0260 = 1.4142135623731*IT_0259;
    const ccomplex_t IT_0261 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0262 = IT_0008*IT_0261;
    const ccomplex_t IT_0263 = 1.4142135623731*IT_0262;
    const ccomplex_t IT_0264 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0265 = IT_0002*IT_0264;
    const ccomplex_t IT_0266 = 1.4142135623731*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*(IT_0260 + -IT_0263 +
       0.333333333333333*IT_0266);
    const ccomplex_t IT_0268 = 3*IT_0267;
    const ccomplex_t IT_0269 = 0.166666666666667*IT_0268;
    const ccomplex_t IT_0270 = IT_0170*IT_0269;
    const ccomplex_t IT_0271 = IT_0168*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = m_b*IT_0272;
    const ccomplex_t IT_0274 = m_b*conj(N_d1)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0275 = IT_0060*IT_0274;
    const ccomplex_t IT_0276 = 1.4142135623731*IT_0275;
    const ccomplex_t IT_0277 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0278 = IT_0002*IT_0277;
    const ccomplex_t IT_0279 = 1.4142135623731*IT_0278;
    const ccomplex_t IT_0280 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0281 = IT_0008*IT_0280;
    const ccomplex_t IT_0282 = 1.4142135623731*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*(IT_0276 + 0.333333333333333
      *IT_0279 + -IT_0282);
    const ccomplex_t IT_0284 = 3*IT_0283;
    const ccomplex_t IT_0285 = 0.166666666666667*IT_0284;
    const ccomplex_t IT_0286 = IT_0069*IT_0285;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*IT_0286;
    const ccomplex_t IT_0288 = m_N_2*IT_0030;
    const ccomplex_t IT_0289 = m_N_1*IT_0052;
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = IT_0287*IT_0290;
    const ccomplex_t IT_0292 = (-2)*U_sb_00*((IT_0000*IT_0008*IT_0059*IT_0061
      *IT_0067 + (-0.5)*m_W*IT_0145*(IT_0008 + 1./3*IT_0007*IT_0146))*conj
      (U_sb_00) + -1./2*IT_0006*IT_0008*IT_0059*IT_0141*conj(U_sb_10)) + U_sb_10
      *(IT_0006*IT_0008*IT_0059*IT_0141*conj(U_sb_00) + (-2)*(IT_0000*IT_0008
      *IT_0059*IT_0061*IT_0067 + -1./3*m_W*IT_0007*IT_0145*IT_0146)*conj(U_sb_10
      ));
    const ccomplex_t IT_0293 = (0 + _Complex_I*-1)*(e_em*IT_0292 + (
      -1.4142135623731)*A_b*IT_0000*IT_0140*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0294 = (-0.5)*IT_0293;
    const ccomplex_t IT_0295 = IT_0269*IT_0294;
    const ccomplex_t IT_0296 = IT_0131*IT_0295;
    const ccomplex_t IT_0297 = (0 + _Complex_I*1)*IT_0296;
    const ccomplex_t IT_0298 = IT_0196*IT_0269;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*IT_0298;
    const ccomplex_t IT_0300 = IT_0195*IT_0299;
    const ccomplex_t IT_0301 = m_b*conj(N_d4)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0302 = IT_0060*IT_0301;
    const ccomplex_t IT_0303 = 1.4142135623731*IT_0302;
    const ccomplex_t IT_0304 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0305 = IT_0002*IT_0304;
    const ccomplex_t IT_0306 = 1.4142135623731*IT_0305;
    const ccomplex_t IT_0307 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0308 = IT_0008*IT_0307;
    const ccomplex_t IT_0309 = 1.4142135623731*IT_0308;
    const ccomplex_t IT_0310 = (0 + _Complex_I*1)*(IT_0303 + 0.333333333333333
      *IT_0306 + -IT_0309);
    const ccomplex_t IT_0311 = 3*IT_0310;
    const ccomplex_t IT_0312 = 0.166666666666667*IT_0311;
    const ccomplex_t IT_0313 = IT_0208*IT_0312;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*IT_0313;
    const ccomplex_t IT_0315 = m_N_2*IT_0231;
    const ccomplex_t IT_0316 = m_N_4*IT_0253;
    const ccomplex_t IT_0317 = IT_0315 + IT_0316;
    const ccomplex_t IT_0318 = IT_0314*IT_0317;
    const ccomplex_t IT_0319 = m_b*conj(N_d3)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0320 = IT_0060*IT_0319;
    const ccomplex_t IT_0321 = 1.4142135623731*IT_0320;
    const ccomplex_t IT_0322 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0323 = IT_0008*IT_0322;
    const ccomplex_t IT_0324 = 1.4142135623731*IT_0323;
    const ccomplex_t IT_0325 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0326 = IT_0002*IT_0325;
    const ccomplex_t IT_0327 = 1.4142135623731*IT_0326;
    const ccomplex_t IT_0328 = (0 + _Complex_I*1)*(IT_0321 + -IT_0324 +
       0.333333333333333*IT_0327);
    const ccomplex_t IT_0329 = 3*IT_0328;
    const ccomplex_t IT_0330 = 0.166666666666667*IT_0329;
    const ccomplex_t IT_0331 = IT_0126*IT_0330;
    const ccomplex_t IT_0332 = (0 + _Complex_I*1)*IT_0331;
    const ccomplex_t IT_0333 = m_N_2*IT_0093;
    const ccomplex_t IT_0334 = m_N_3*IT_0115;
    const ccomplex_t IT_0335 = IT_0333 + IT_0334;
    const ccomplex_t IT_0336 = IT_0332*IT_0335;
    const ccomplex_t IT_0337 = m_b*conj(N_d2)*e_em*IT_0061*U_sb_11;
    const ccomplex_t IT_0338 = IT_0060*IT_0337;
    const ccomplex_t IT_0339 = 1.4142135623731*IT_0338;
    const ccomplex_t IT_0340 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0341 = IT_0008*IT_0340;
    const ccomplex_t IT_0342 = 1.4142135623731*IT_0341;
    const ccomplex_t IT_0343 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0344 = IT_0002*IT_0343;
    const ccomplex_t IT_0345 = 1.4142135623731*IT_0344;
    const ccomplex_t IT_0346 = (0 + _Complex_I*1)*(IT_0339 + -IT_0342 +
       0.333333333333333*IT_0345);
    const ccomplex_t IT_0347 = 3*IT_0346;
    const ccomplex_t IT_0348 = 0.166666666666667*IT_0347;
    const ccomplex_t IT_0349 = (-2)*U_sb_00*((IT_0000*IT_0008*IT_0059*IT_0061
      *IT_0067 + (-0.5)*m_W*IT_0145*(IT_0008 + 1./3*IT_0007*IT_0146))*conj
      (U_sb_01) + -1./2*IT_0006*IT_0008*IT_0059*IT_0141*conj(U_sb_11)) + U_sb_10
      *(IT_0006*IT_0008*IT_0059*IT_0141*conj(U_sb_01) + (-2)*(IT_0000*IT_0008
      *IT_0059*IT_0061*IT_0067 + -1./3*m_W*IT_0007*IT_0145*IT_0146)*conj(U_sb_11
      ));
    const ccomplex_t IT_0350 = (0 + _Complex_I*-1)*(e_em*IT_0349 + (
      -1.4142135623731)*A_b*IT_0000*IT_0140*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0351 = (-0.5)*IT_0350;
    const ccomplex_t IT_0352 = IT_0348*IT_0351;
    const ccomplex_t IT_0353 = IT_0153*IT_0352;
    const ccomplex_t IT_0354 = (0 + _Complex_I*1)*IT_0353;
    const ccomplex_t IT_0355 = 2*IT_0273 + -IT_0291 + -IT_0297 + -IT_0300 + 
      -IT_0318 + -IT_0336 + -IT_0354;
    const ccomplex_t IT_0356 = 6*s_13;
    const ccomplex_t IT_0357 = m_b*m_N_2;
    const ccomplex_t IT_0358 = 6*IT_0257;
    const ccomplex_t IT_0359 = IT_0231*IT_0314;
    const ccomplex_t IT_0360 = IT_0183*IT_0299;
    const ccomplex_t IT_0361 = IT_0030*IT_0287;
    const ccomplex_t IT_0362 = IT_0093*IT_0332;
    const ccomplex_t IT_0363 = -IT_0173 + IT_0359 + IT_0360 + IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = m_b*s_14;
    const ccomplex_t IT_0365 = IT_0193*IT_0198;
    const ccomplex_t IT_0366 = IT_0115*IT_0128;
    const ccomplex_t IT_0367 = IT_0210*IT_0253;
    const ccomplex_t IT_0368 = IT_0052*IT_0071;
    const ccomplex_t IT_0369 = IT_0272 + -IT_0365 + -IT_0366 + -IT_0367 + 
      -IT_0368;
    const ccomplex_t IT_0370 = s_34*m_N_2;
    const ccomplex_t IT_0371 = 6*conj(IT_0257);
    const ccomplex_t IT_0372 = 6*conj(IT_0355);
    const ccomplex_t IT_0373 = pow(m_H0, 2);
    const ccomplex_t IT_0374 = m_b*m_N_2*IT_0373;
    const ccomplex_t IT_0375 = s_13*IT_0373;
    const ccomplex_t IT_0376 = s_14*s_34;
    const ccomplex_t IT_0377 = (-2)*IT_0376;
    const ccomplex_t IT_0378 = IT_0375 + IT_0377;
    const ccomplex_t IT_0379 = 6*IT_0355;
    const ccomplex_t IT_0380 = 0.166666666666667*IT_0370;
    const ccomplex_t IT_0381 = (IT_0257*conj(IT_0257) + IT_0355*conj(IT_0355))
      *IT_0356 + conj(IT_0355)*IT_0357*IT_0358 + (IT_0355*IT_0357 + IT_0363
      *IT_0364 + IT_0369*IT_0370)*IT_0371 + (IT_0364*IT_0369 + IT_0363*IT_0370)
      *IT_0372 + 6*conj(IT_0369)*(IT_0363*IT_0374 + -IT_0369*IT_0378 +
       0.166666666666667*IT_0364*IT_0379 + IT_0358*IT_0380) + 6*conj(IT_0363)*
      (0.166666666666667*IT_0358*IT_0364 + IT_0369*IT_0374 + -IT_0363*IT_0378 +
       IT_0379*IT_0380);
    return create_ccomplex_return(IT_0381);
}

