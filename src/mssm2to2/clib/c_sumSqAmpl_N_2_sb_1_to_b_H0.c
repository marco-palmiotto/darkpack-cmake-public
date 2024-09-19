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
    const ccomplex_t IT_0004 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0013 = IT_0003*IT_0012;
    const ccomplex_t IT_0014 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = IT_0005 + IT_0011 + IT_0013 + IT_0015;
    const ccomplex_t IT_0017 = IT_0002*IT_0006;
    const ccomplex_t IT_0018 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0000*IT_0008;
    const ccomplex_t IT_0021 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0022 + -IT_0024 + -IT_0026;
    const ccomplex_t IT_0028 = IT_0016 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = m_N_3*IT_0030;
    const ccomplex_t IT_0032 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0035 = IT_0009*IT_0034;
    const ccomplex_t IT_0036 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0039 = IT_0009*IT_0038;
    const ccomplex_t IT_0040 = IT_0033 + IT_0035 + IT_0037 + IT_0039;
    const ccomplex_t IT_0041 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0042 = IT_0017*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0044 = IT_0020*IT_0043;
    const ccomplex_t IT_0045 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0046 = IT_0017*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0048 = IT_0020*IT_0047;
    const ccomplex_t IT_0049 = -IT_0042 + -IT_0044 + -IT_0046 + -IT_0048;
    const ccomplex_t IT_0050 = IT_0040 + IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = m_N_2*IT_0052;
    const ccomplex_t IT_0054 = IT_0031 + IT_0053;
    const ccomplex_t IT_0055 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0056 = IT_0002*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = cos(beta);
    const ccomplex_t IT_0059 = cpow(IT_0058, -1);
    const ccomplex_t IT_0060 = IT_0008*IT_0059;
    const ccomplex_t IT_0061 = pow(m_W, -1);
    const ccomplex_t IT_0062 = m_b*N_d3*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0057 + 1.5*IT_0064);
    const ccomplex_t IT_0066 = 0.666666666666667*IT_0065;
    const ccomplex_t IT_0067 = 1.5*IT_0066;
    const ccomplex_t IT_0068 = (-0.333333333333333)*IT_0067;
    const ccomplex_t IT_0069 = pow(m_b, 2);
    const ccomplex_t IT_0070 = pow(m_sb_1, 2);
    const ccomplex_t IT_0071 = cpow((-2)*s_23 + IT_0069 + IT_0070 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0068*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0054*IT_0073;
    const ccomplex_t IT_0075 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0076 = IT_0002*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = m_b*N_d4*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0079 = IT_0060*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 1.5*IT_0080);
    const ccomplex_t IT_0082 = (-0.333333333333333)*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_23 + IT_0069 + IT_0070 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0087 = IT_0003*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0089 = IT_0009*IT_0088;
    const ccomplex_t IT_0090 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0091 = IT_0003*IT_0090;
    const ccomplex_t IT_0092 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0093 = IT_0009*IT_0092;
    const ccomplex_t IT_0094 = IT_0087 + IT_0089 + IT_0091 + IT_0093;
    const ccomplex_t IT_0095 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0096 = IT_0017*IT_0095;
    const ccomplex_t IT_0097 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0098 = IT_0020*IT_0097;
    const ccomplex_t IT_0099 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0100 = IT_0017*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0102 = IT_0020*IT_0101;
    const ccomplex_t IT_0103 = -IT_0096 + -IT_0098 + -IT_0100 + -IT_0102;
    const ccomplex_t IT_0104 = IT_0094 + IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = m_N_2*IT_0106;
    const ccomplex_t IT_0108 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0109 = IT_0003*IT_0108;
    const ccomplex_t IT_0110 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0111 = IT_0009*IT_0110;
    const ccomplex_t IT_0112 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0113 = IT_0003*IT_0112;
    const ccomplex_t IT_0114 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0115 = IT_0009*IT_0114;
    const ccomplex_t IT_0116 = IT_0109 + IT_0111 + IT_0113 + IT_0115;
    const ccomplex_t IT_0117 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0118 = IT_0017*IT_0117;
    const ccomplex_t IT_0119 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0120 = IT_0020*IT_0119;
    const ccomplex_t IT_0121 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0122 = IT_0017*IT_0121;
    const ccomplex_t IT_0123 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0124 = IT_0020*IT_0123;
    const ccomplex_t IT_0125 = -IT_0118 + -IT_0120 + -IT_0122 + -IT_0124;
    const ccomplex_t IT_0126 = IT_0116 + IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = m_N_4*IT_0128;
    const ccomplex_t IT_0130 = IT_0107 + IT_0129;
    const ccomplex_t IT_0131 = IT_0085*IT_0130;
    const ccomplex_t IT_0132 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0133 = IT_0003*IT_0132;
    const ccomplex_t IT_0134 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0135 = IT_0009*IT_0134;
    const ccomplex_t IT_0136 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0137 = IT_0003*IT_0136;
    const ccomplex_t IT_0138 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0139 = IT_0009*IT_0138;
    const ccomplex_t IT_0140 = IT_0133 + IT_0135 + IT_0137 + IT_0139;
    const ccomplex_t IT_0141 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0142 = IT_0017*IT_0141;
    const ccomplex_t IT_0143 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0144 = IT_0020*IT_0143;
    const ccomplex_t IT_0145 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0146 = IT_0017*IT_0145;
    const ccomplex_t IT_0147 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0148 = IT_0020*IT_0147;
    const ccomplex_t IT_0149 = -IT_0142 + -IT_0144 + -IT_0146 + -IT_0148;
    const ccomplex_t IT_0150 = IT_0140 + IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = m_N_1*IT_0152;
    const ccomplex_t IT_0154 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0155 = IT_0003*IT_0154;
    const ccomplex_t IT_0156 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0157 = IT_0009*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0159 = IT_0003*IT_0158;
    const ccomplex_t IT_0160 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0161 = IT_0009*IT_0160;
    const ccomplex_t IT_0162 = IT_0155 + IT_0157 + IT_0159 + IT_0161;
    const ccomplex_t IT_0163 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0164 = IT_0017*IT_0163;
    const ccomplex_t IT_0165 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0166 = IT_0020*IT_0165;
    const ccomplex_t IT_0167 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0168 = IT_0017*IT_0167;
    const ccomplex_t IT_0169 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0170 = IT_0020*IT_0169;
    const ccomplex_t IT_0171 = -IT_0164 + -IT_0166 + -IT_0168 + -IT_0170;
    const ccomplex_t IT_0172 = IT_0162 + IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0173;
    const ccomplex_t IT_0175 = m_N_2*IT_0174;
    const ccomplex_t IT_0176 = IT_0153 + IT_0175;
    const ccomplex_t IT_0177 = m_b*N_d1*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0178 = IT_0060*IT_0177;
    const ccomplex_t IT_0179 = 1.4142135623731*IT_0178;
    const ccomplex_t IT_0180 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0181 = IT_0002*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0179 + 0.666666666666667
      *IT_0182);
    const ccomplex_t IT_0184 = 1.5*IT_0183;
    const ccomplex_t IT_0185 = (-0.333333333333333)*IT_0184;
    const ccomplex_t IT_0186 = cpow((-2)*s_23 + IT_0069 + IT_0070 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0187 = IT_0185*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = IT_0176*IT_0188;
    const ccomplex_t IT_0190 = pow(m_N_2, 2);
    const ccomplex_t IT_0191 = cpow((-2)*s_13 + IT_0069 + IT_0190 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0192 = cpow(V_tb, 2);
    const ccomplex_t IT_0193 = m_b*mu_h*IT_0061;
    const ccomplex_t IT_0194 = sin(beta);
    const ccomplex_t IT_0195 = -IT_0006*IT_0194;
    const ccomplex_t IT_0196 = IT_0000*IT_0058;
    const ccomplex_t IT_0197 = IT_0195 + IT_0196;
    const ccomplex_t IT_0198 = cpow(IT_0001, -2);
    const ccomplex_t IT_0199 = (-2)*U_sb_00*((IT_0000*IT_0008*IT_0059*IT_0061
      *IT_0069 + (-0.5)*m_W*IT_0197*(IT_0008 + 1./3*IT_0007*IT_0198))*conj
      (U_sb_00) + -1./2*IT_0006*IT_0008*IT_0059*IT_0193*conj(U_sb_10)) + U_sb_10
      *(IT_0006*IT_0008*IT_0059*IT_0193*conj(U_sb_00) + (-2)*(IT_0000*IT_0008
      *IT_0059*IT_0061*IT_0069 + -1./3*m_W*IT_0007*IT_0197*IT_0198)*conj(U_sb_10
      ));
    const ccomplex_t IT_0200 = (0 + _Complex_I*-1)*(e_em*IT_0199 + (
      -1.4142135623731)*A_b*IT_0000*IT_0192*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0201 = (-0.5)*IT_0200;
    const ccomplex_t IT_0202 = m_b*N_d2*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0203 = IT_0060*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0206 = IT_0002*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0204 + 0.666666666666667
      *IT_0207);
    const ccomplex_t IT_0209 = 1.5*IT_0208;
    const ccomplex_t IT_0210 = 0.666666666666667*IT_0209;
    const ccomplex_t IT_0211 = 1.5*IT_0210;
    const ccomplex_t IT_0212 = (-0.333333333333333)*IT_0211;
    const ccomplex_t IT_0213 = IT_0201*IT_0212;
    const ccomplex_t IT_0214 = IT_0191*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = cpow((-2)*s_13 + IT_0069 + IT_0190 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0217 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0218 = IT_0002*IT_0217;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = m_b*N_d2*e_em*IT_0061*U_sb_01;
    const ccomplex_t IT_0221 = IT_0060*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*(IT_0219 + 1.5*IT_0222);
    const ccomplex_t IT_0224 = (-0.333333333333333)*IT_0223;
    const ccomplex_t IT_0225 = (-2)*U_sb_00*((IT_0000*IT_0008*IT_0059*IT_0061
      *IT_0069 + (-0.5)*m_W*IT_0197*(IT_0008 + 1./3*IT_0007*IT_0198))*conj
      (U_sb_01) + -1./2*IT_0006*IT_0008*IT_0059*IT_0193*conj(U_sb_11)) + U_sb_10
      *(IT_0006*IT_0008*IT_0059*IT_0193*conj(U_sb_01) + (-2)*(IT_0000*IT_0008
      *IT_0059*IT_0061*IT_0069 + -1./3*m_W*IT_0007*IT_0197*IT_0198)*conj(U_sb_11
      ));
    const ccomplex_t IT_0226 = (0 + _Complex_I*-1)*(e_em*IT_0225 + (
      -1.4142135623731)*A_b*IT_0000*IT_0192*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0227 = (-0.5)*IT_0226;
    const ccomplex_t IT_0228 = IT_0224*IT_0227;
    const ccomplex_t IT_0229 = IT_0216*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = cpow((-2)*s_12 + IT_0069 + -IT_0070 + -IT_0190 
      + -reg_prop, -1);
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0008
      *IT_0059*IT_0061;
    const ccomplex_t IT_0233 = (-0.5)*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*(1.5*IT_0204 + IT_0207);
    const ccomplex_t IT_0235 = 0.666666666666667*IT_0234;
    const ccomplex_t IT_0236 = 1.5*IT_0235;
    const ccomplex_t IT_0237 = (-0.333333333333333)*IT_0236;
    const ccomplex_t IT_0238 = IT_0233*IT_0237;
    const ccomplex_t IT_0239 = IT_0231*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*IT_0239;
    const ccomplex_t IT_0241 = m_b*IT_0240;
    const ccomplex_t IT_0242 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0243 = IT_0003*IT_0242;
    const ccomplex_t IT_0244 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0245 = IT_0017*IT_0244;
    const ccomplex_t IT_0246 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0247 = IT_0020*IT_0246;
    const ccomplex_t IT_0248 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0249 = IT_0009*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*(IT_0243 + -IT_0245 + 
      -IT_0247 + IT_0249);
    const ccomplex_t IT_0251 = m_N_2*IT_0250;
    const ccomplex_t IT_0252 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0253 = IT_0003*IT_0252;
    const ccomplex_t IT_0254 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0255 = IT_0017*IT_0254;
    const ccomplex_t IT_0256 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0257 = IT_0020*IT_0256;
    const ccomplex_t IT_0258 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0259 = IT_0009*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*(IT_0253 + -IT_0255 + 
      -IT_0257 + IT_0259);
    const ccomplex_t IT_0261 = m_N_2*IT_0260;
    const ccomplex_t IT_0262 = IT_0251 + IT_0261;
    const ccomplex_t IT_0263 = cpow((-2)*s_23 + IT_0069 + IT_0070 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0264 = IT_0212*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = IT_0262*IT_0265;
    const ccomplex_t IT_0267 = IT_0074 + IT_0131 + IT_0189 + IT_0215 + IT_0230
       + (-2)*IT_0241 + IT_0266;
    const ccomplex_t IT_0268 = (-2)*U_sb_00*((IT_0000*IT_0008*IT_0059*IT_0061
      *IT_0069 + (-0.5)*m_W*IT_0197*(IT_0008 + 1./3*IT_0007*IT_0198))*conj
      (U_sb_00) + -1./2*IT_0006*IT_0008*IT_0059*IT_0193*conj(U_sb_10)) + U_sb_10
      *(IT_0006*IT_0008*IT_0059*IT_0193*conj(U_sb_00) + (-2)*(IT_0000*IT_0008
      *IT_0059*IT_0061*IT_0069 + -1./3*m_W*IT_0007*IT_0197*IT_0198)*conj(U_sb_10
      ));
    const ccomplex_t IT_0269 = (0 + _Complex_I*-1)*(e_em*IT_0268 + (
      -1.4142135623731)*A_b*IT_0000*IT_0192*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0270 = (-0.5)*IT_0269;
    const ccomplex_t IT_0271 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0272 = IT_0002*IT_0271;
    const ccomplex_t IT_0273 = 1.4142135623731*IT_0272;
    const ccomplex_t IT_0274 = m_b*conj(N_d2)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0275 = IT_0060*IT_0274;
    const ccomplex_t IT_0276 = 1.4142135623731*IT_0275;
    const ccomplex_t IT_0277 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0278 = IT_0008*IT_0277;
    const ccomplex_t IT_0279 = 1.4142135623731*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*(IT_0273 + 3*IT_0276 + (-3)
      *IT_0279);
    const ccomplex_t IT_0281 = 0.333333333333333*IT_0280;
    const ccomplex_t IT_0282 = 3*IT_0281;
    const ccomplex_t IT_0283 = (-0.333333333333333)*IT_0282;
    const ccomplex_t IT_0284 = -IT_0283;
    const ccomplex_t IT_0285 = 3*IT_0284;
    const ccomplex_t IT_0286 = 0.166666666666667*IT_0285;
    const ccomplex_t IT_0287 = IT_0270*IT_0286;
    const ccomplex_t IT_0288 = IT_0191*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0291 = IT_0008*IT_0290;
    const ccomplex_t IT_0292 = 1.4142135623731*IT_0291;
    const ccomplex_t IT_0293 = m_b*conj(N_d2)*e_em*IT_0061*U_sb_11;
    const ccomplex_t IT_0294 = IT_0060*IT_0293;
    const ccomplex_t IT_0295 = 1.4142135623731*IT_0294;
    const ccomplex_t IT_0296 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0297 = IT_0002*IT_0296;
    const ccomplex_t IT_0298 = 1.4142135623731*IT_0297;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*(IT_0292 + -IT_0295 + (
      -0.333333333333333)*IT_0298);
    const ccomplex_t IT_0300 = (-3)*IT_0299;
    const ccomplex_t IT_0301 = 0.166666666666667*IT_0300;
    const ccomplex_t IT_0302 = (-2)*U_sb_00*((IT_0000*IT_0008*IT_0059*IT_0061
      *IT_0069 + (-0.5)*m_W*IT_0197*(IT_0008 + 1./3*IT_0007*IT_0198))*conj
      (U_sb_01) + -1./2*IT_0006*IT_0008*IT_0059*IT_0193*conj(U_sb_11)) + U_sb_10
      *(IT_0006*IT_0008*IT_0059*IT_0193*conj(U_sb_01) + (-2)*(IT_0000*IT_0008
      *IT_0059*IT_0061*IT_0069 + -1./3*m_W*IT_0007*IT_0197*IT_0198)*conj(U_sb_11
      ));
    const ccomplex_t IT_0303 = (0 + _Complex_I*-1)*(e_em*IT_0302 + (
      -1.4142135623731)*A_b*IT_0000*IT_0192*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0304 = (-0.5)*IT_0303;
    const ccomplex_t IT_0305 = IT_0301*IT_0304;
    const ccomplex_t IT_0306 = IT_0216*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*IT_0306;
    const ccomplex_t IT_0308 = IT_0233*IT_0286;
    const ccomplex_t IT_0309 = IT_0231*IT_0308;
    const ccomplex_t IT_0310 = (0 + _Complex_I*1)*IT_0309;
    const ccomplex_t IT_0311 = m_b*IT_0310;
    const ccomplex_t IT_0312 = IT_0263*IT_0286;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*IT_0312;
    const ccomplex_t IT_0314 = IT_0262*IT_0313;
    const ccomplex_t IT_0315 = m_N_2*IT_0152;
    const ccomplex_t IT_0316 = m_N_1*IT_0174;
    const ccomplex_t IT_0317 = IT_0315 + IT_0316;
    const ccomplex_t IT_0318 = m_b*conj(N_d1)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0319 = IT_0060*IT_0318;
    const ccomplex_t IT_0320 = 1.4142135623731*IT_0319;
    const ccomplex_t IT_0321 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0322 = IT_0008*IT_0321;
    const ccomplex_t IT_0323 = 1.4142135623731*IT_0322;
    const ccomplex_t IT_0324 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0325 = IT_0002*IT_0324;
    const ccomplex_t IT_0326 = 1.4142135623731*IT_0325;
    const ccomplex_t IT_0327 = (0 + _Complex_I*1)*(IT_0320 + -IT_0323 +
       0.333333333333333*IT_0326);
    const ccomplex_t IT_0328 = -IT_0327;
    const ccomplex_t IT_0329 = (-3)*IT_0328;
    const ccomplex_t IT_0330 = (-0.333333333333333)*IT_0329;
    const ccomplex_t IT_0331 = -IT_0330;
    const ccomplex_t IT_0332 = 3*IT_0331;
    const ccomplex_t IT_0333 = 0.166666666666667*IT_0332;
    const ccomplex_t IT_0334 = IT_0186*IT_0333;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = IT_0317*IT_0335;
    const ccomplex_t IT_0337 = m_N_2*IT_0030;
    const ccomplex_t IT_0338 = m_N_3*IT_0052;
    const ccomplex_t IT_0339 = IT_0337 + IT_0338;
    const ccomplex_t IT_0340 = m_b*conj(N_d3)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0341 = IT_0060*IT_0340;
    const ccomplex_t IT_0342 = 1.4142135623731*IT_0341;
    const ccomplex_t IT_0343 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0344 = IT_0008*IT_0343;
    const ccomplex_t IT_0345 = 1.4142135623731*IT_0344;
    const ccomplex_t IT_0346 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0347 = IT_0002*IT_0346;
    const ccomplex_t IT_0348 = 1.4142135623731*IT_0347;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*(IT_0342 + -IT_0345 +
       0.333333333333333*IT_0348);
    const ccomplex_t IT_0350 = 3*IT_0349;
    const ccomplex_t IT_0351 = 0.166666666666667*IT_0350;
    const ccomplex_t IT_0352 = IT_0071*IT_0351;
    const ccomplex_t IT_0353 = (0 + _Complex_I*1)*IT_0352;
    const ccomplex_t IT_0354 = IT_0339*IT_0353;
    const ccomplex_t IT_0355 = m_b*conj(N_d4)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0356 = IT_0060*IT_0355;
    const ccomplex_t IT_0357 = 1.4142135623731*IT_0356;
    const ccomplex_t IT_0358 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0359 = IT_0008*IT_0358;
    const ccomplex_t IT_0360 = 1.4142135623731*IT_0359;
    const ccomplex_t IT_0361 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0362 = IT_0002*IT_0361;
    const ccomplex_t IT_0363 = 1.4142135623731*IT_0362;
    const ccomplex_t IT_0364 = (0 + _Complex_I*1)*(IT_0357 + -IT_0360 +
       0.333333333333333*IT_0363);
    const ccomplex_t IT_0365 = 3*IT_0364;
    const ccomplex_t IT_0366 = 0.166666666666667*IT_0365;
    const ccomplex_t IT_0367 = IT_0083*IT_0366;
    const ccomplex_t IT_0368 = (0 + _Complex_I*1)*IT_0367;
    const ccomplex_t IT_0369 = m_N_2*IT_0128;
    const ccomplex_t IT_0370 = m_N_4*IT_0106;
    const ccomplex_t IT_0371 = IT_0369 + IT_0370;
    const ccomplex_t IT_0372 = IT_0368*IT_0371;
    const ccomplex_t IT_0373 = -IT_0289 + -IT_0307 + 2*IT_0311 + -IT_0314 + 
      -IT_0336 + -IT_0354 + -IT_0372;
    const ccomplex_t IT_0374 = 6*s_13;
    const ccomplex_t IT_0375 = m_b*m_N_2;
    const ccomplex_t IT_0376 = 6*IT_0267;
    const ccomplex_t IT_0377 = IT_0128*IT_0368;
    const ccomplex_t IT_0378 = IT_0152*IT_0335;
    const ccomplex_t IT_0379 = IT_0030*IT_0353;
    const ccomplex_t IT_0380 = IT_0250*IT_0313;
    const ccomplex_t IT_0381 = -IT_0240 + IT_0377 + IT_0378 + IT_0379 + IT_0380;
    const ccomplex_t IT_0382 = m_b*s_14;
    const ccomplex_t IT_0383 = IT_0237*IT_0263;
    const ccomplex_t IT_0384 = (0 + _Complex_I*1)*IT_0383;
    const ccomplex_t IT_0385 = IT_0260*IT_0384;
    const ccomplex_t IT_0386 = (0 + _Complex_I*1)*(1.5*IT_0179 + IT_0182);
    const ccomplex_t IT_0387 = (-0.333333333333333)*IT_0386;
    const ccomplex_t IT_0388 = IT_0186*IT_0387;
    const ccomplex_t IT_0389 = (0 + _Complex_I*1)*IT_0388;
    const ccomplex_t IT_0390 = IT_0174*IT_0389;
    const ccomplex_t IT_0391 = IT_0052*IT_0073;
    const ccomplex_t IT_0392 = IT_0085*IT_0106;
    const ccomplex_t IT_0393 = IT_0310 + -IT_0385 + -IT_0390 + -IT_0391 + 
      -IT_0392;
    const ccomplex_t IT_0394 = s_34*m_N_2;
    const ccomplex_t IT_0395 = 6*conj(IT_0267);
    const ccomplex_t IT_0396 = 6*conj(IT_0373);
    const ccomplex_t IT_0397 = pow(m_H0, 2);
    const ccomplex_t IT_0398 = m_b*m_N_2*IT_0397;
    const ccomplex_t IT_0399 = s_13*IT_0397;
    const ccomplex_t IT_0400 = s_14*s_34;
    const ccomplex_t IT_0401 = (-2)*IT_0400;
    const ccomplex_t IT_0402 = IT_0399 + IT_0401;
    const ccomplex_t IT_0403 = 6*IT_0373;
    const ccomplex_t IT_0404 = 0.166666666666667*IT_0394;
    const ccomplex_t IT_0405 = (IT_0267*conj(IT_0267) + IT_0373*conj(IT_0373))
      *IT_0374 + conj(IT_0373)*IT_0375*IT_0376 + (IT_0373*IT_0375 + IT_0381
      *IT_0382 + IT_0393*IT_0394)*IT_0395 + (IT_0382*IT_0393 + IT_0381*IT_0394)
      *IT_0396 + 6*conj(IT_0393)*(IT_0381*IT_0398 + -IT_0393*IT_0402 +
       0.166666666666667*IT_0382*IT_0403 + IT_0376*IT_0404) + 6*conj(IT_0381)*
      (0.166666666666667*IT_0376*IT_0382 + IT_0393*IT_0398 + -IT_0381*IT_0402 +
       IT_0403*IT_0404);
    return create_ccomplex_return(IT_0405);
}

