#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_sb_1_to_b_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_sb_1_to_b_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
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
    const ccomplex_t IT_0004 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = IT_0005 + IT_0009 + IT_0011 + IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0000*IT_0016;
    const ccomplex_t IT_0018 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0006*IT_0016;
    const ccomplex_t IT_0021 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0022 + -IT_0024 + -IT_0026;
    const ccomplex_t IT_0028 = IT_0014 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = m_N_2*IT_0030;
    const ccomplex_t IT_0032 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0035 = IT_0007*IT_0034;
    const ccomplex_t IT_0036 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0039 = IT_0007*IT_0038;
    const ccomplex_t IT_0040 = IT_0033 + IT_0035 + IT_0037 + IT_0039;
    const ccomplex_t IT_0041 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0042 = IT_0017*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0044 = IT_0020*IT_0043;
    const ccomplex_t IT_0045 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0046 = IT_0017*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0048 = IT_0020*IT_0047;
    const ccomplex_t IT_0049 = -IT_0042 + -IT_0044 + -IT_0046 + -IT_0048;
    const ccomplex_t IT_0050 = IT_0040 + IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = m_N_3*IT_0052;
    const ccomplex_t IT_0054 = IT_0031 + IT_0053;
    const ccomplex_t IT_0055 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0056 = IT_0002*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = cos(beta);
    const ccomplex_t IT_0059 = cpow(IT_0058, -1);
    const ccomplex_t IT_0060 = IT_0016*IT_0059;
    const ccomplex_t IT_0061 = pow(m_W, -1);
    const ccomplex_t IT_0062 = m_b*N_d2*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0057 + 1.5*IT_0064);
    const ccomplex_t IT_0066 = 0.666666666666667*IT_0065;
    const ccomplex_t IT_0067 = 1.5*IT_0066;
    const ccomplex_t IT_0068 = 0.666666666666667*IT_0067;
    const ccomplex_t IT_0069 = 1.5*IT_0068;
    const ccomplex_t IT_0070 = 0.666666666666667*IT_0069;
    const ccomplex_t IT_0071 = 1.5*IT_0070;
    const ccomplex_t IT_0072 = 0.666666666666667*IT_0071;
    const ccomplex_t IT_0073 = 1.5*IT_0072;
    const ccomplex_t IT_0074 = 0.666666666666667*IT_0073;
    const ccomplex_t IT_0075 = 1.5*IT_0074;
    const ccomplex_t IT_0076 = 0.666666666666667*IT_0075;
    const ccomplex_t IT_0077 = 1.5*IT_0076;
    const ccomplex_t IT_0078 = (-0.333333333333333)*IT_0077;
    const ccomplex_t IT_0079 = pow(m_b, 2);
    const ccomplex_t IT_0080 = pow(m_sb_1, 2);
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0079 + IT_0080 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0078*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = IT_0054*IT_0083;
    const ccomplex_t IT_0085 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0086 = IT_0003*IT_0085;
    const ccomplex_t IT_0087 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0088 = IT_0007*IT_0087;
    const ccomplex_t IT_0089 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0090 = IT_0003*IT_0089;
    const ccomplex_t IT_0091 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0092 = IT_0007*IT_0091;
    const ccomplex_t IT_0093 = IT_0086 + IT_0088 + IT_0090 + IT_0092;
    const ccomplex_t IT_0094 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0095 = IT_0017*IT_0094;
    const ccomplex_t IT_0096 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0097 = IT_0020*IT_0096;
    const ccomplex_t IT_0098 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0099 = IT_0017*IT_0098;
    const ccomplex_t IT_0100 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0101 = IT_0020*IT_0100;
    const ccomplex_t IT_0102 = -IT_0095 + -IT_0097 + -IT_0099 + -IT_0101;
    const ccomplex_t IT_0103 = IT_0093 + IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = m_N_4*IT_0105;
    const ccomplex_t IT_0107 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0108 = IT_0003*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0110 = IT_0007*IT_0109;
    const ccomplex_t IT_0111 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0112 = IT_0003*IT_0111;
    const ccomplex_t IT_0113 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0114 = IT_0007*IT_0113;
    const ccomplex_t IT_0115 = IT_0108 + IT_0110 + IT_0112 + IT_0114;
    const ccomplex_t IT_0116 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0117 = IT_0017*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0119 = IT_0020*IT_0118;
    const ccomplex_t IT_0120 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0121 = IT_0017*IT_0120;
    const ccomplex_t IT_0122 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0123 = IT_0020*IT_0122;
    const ccomplex_t IT_0124 = -IT_0117 + -IT_0119 + -IT_0121 + -IT_0123;
    const ccomplex_t IT_0125 = IT_0115 + IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = (-0.5)*IT_0126;
    const ccomplex_t IT_0128 = m_N_3*IT_0127;
    const ccomplex_t IT_0129 = IT_0106 + IT_0128;
    const ccomplex_t IT_0130 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0131 = IT_0002*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = m_b*N_d4*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0134 = IT_0060*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + 1.5*IT_0135);
    const ccomplex_t IT_0137 = (-0.333333333333333)*IT_0136;
    const ccomplex_t IT_0138 = cpow((-2)*s_23 + IT_0079 + IT_0080 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = IT_0129*IT_0140;
    const ccomplex_t IT_0142 = pow(m_N_3, 2);
    const ccomplex_t IT_0143 = cpow((-2)*s_13 + IT_0079 + IT_0142 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0144 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0145 = IT_0002*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = m_b*N_d3*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0148 = IT_0060*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0146 + 1.5*IT_0149);
    const ccomplex_t IT_0151 = 0.666666666666667*IT_0150;
    const ccomplex_t IT_0152 = 1.5*IT_0151;
    const ccomplex_t IT_0153 = (-0.333333333333333)*IT_0152;
    const ccomplex_t IT_0154 = cpow(V_tb, 2);
    const ccomplex_t IT_0155 = m_b*mu_h*IT_0061;
    const ccomplex_t IT_0156 = IT_0006*IT_0058;
    const ccomplex_t IT_0157 = sin(beta);
    const ccomplex_t IT_0158 = IT_0000*IT_0157;
    const ccomplex_t IT_0159 = IT_0156 + IT_0158;
    const ccomplex_t IT_0160 = cpow(IT_0001, -2);
    const ccomplex_t IT_0161 = 2*U_sb_00*((IT_0006*IT_0016*IT_0059*IT_0061
      *IT_0079 + (-0.5)*m_W*IT_0159*(IT_0016 + 1./3*IT_0015*IT_0160))*conj
      (U_sb_00) + 1./2*IT_0000*IT_0016*IT_0059*IT_0155*conj(U_sb_10)) + U_sb_10*
      (IT_0000*IT_0016*IT_0059*IT_0155*conj(U_sb_00) + 2*(IT_0006*IT_0016
      *IT_0059*IT_0061*IT_0079 + -1./3*m_W*IT_0015*IT_0159*IT_0160)*conj(U_sb_10
      ));
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*(e_em*IT_0161 +
       1.4142135623731*A_b*IT_0006*IT_0154*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0163 = 0.5*IT_0162;
    const ccomplex_t IT_0164 = IT_0153*IT_0163;
    const ccomplex_t IT_0165 = IT_0143*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = cpow((-2)*s_13 + IT_0079 + IT_0142 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0168 = m_b*N_d3*e_em*IT_0061*U_sb_01;
    const ccomplex_t IT_0169 = IT_0060*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0172 = IT_0002*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0170 + 0.666666666666667
      *IT_0173);
    const ccomplex_t IT_0175 = 1.5*IT_0174;
    const ccomplex_t IT_0176 = 0.666666666666667*IT_0175;
    const ccomplex_t IT_0177 = 1.5*IT_0176;
    const ccomplex_t IT_0178 = (-0.333333333333333)*IT_0177;
    const ccomplex_t IT_0179 = 2*U_sb_00*((IT_0006*IT_0016*IT_0059*IT_0061
      *IT_0079 + (-0.5)*m_W*IT_0159*(IT_0016 + 1./3*IT_0015*IT_0160))*conj
      (U_sb_01) + 1./2*IT_0000*IT_0016*IT_0059*IT_0155*conj(U_sb_11)) + U_sb_10*
      (IT_0000*IT_0016*IT_0059*IT_0155*conj(U_sb_01) + 2*(IT_0006*IT_0016
      *IT_0059*IT_0061*IT_0079 + -1./3*m_W*IT_0015*IT_0159*IT_0160)*conj(U_sb_11
      ));
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*(e_em*IT_0179 +
       1.4142135623731*A_b*IT_0006*IT_0154*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0181 = 0.5*IT_0180;
    const ccomplex_t IT_0182 = IT_0178*IT_0181;
    const ccomplex_t IT_0183 = IT_0167*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = cpow((-2)*s_12 + IT_0079 + -IT_0080 + -IT_0142 
      + -reg_prop, -1);
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*m_b*e_em*IT_0006*IT_0016
      *IT_0059*IT_0061;
    const ccomplex_t IT_0187 = 0.5*IT_0186;
    const ccomplex_t IT_0188 = IT_0153*IT_0187;
    const ccomplex_t IT_0189 = IT_0185*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = m_b*IT_0190;
    const ccomplex_t IT_0192 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0193 = IT_0003*IT_0192;
    const ccomplex_t IT_0194 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0195 = IT_0007*IT_0194;
    const ccomplex_t IT_0196 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0197 = IT_0017*IT_0196;
    const ccomplex_t IT_0198 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0199 = IT_0020*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*(IT_0193 + IT_0195 + 
      -IT_0197 + -IT_0199);
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = m_N_3*IT_0201;
    const ccomplex_t IT_0203 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0204 = IT_0003*IT_0203;
    const ccomplex_t IT_0205 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0206 = IT_0007*IT_0205;
    const ccomplex_t IT_0207 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0208 = IT_0017*IT_0207;
    const ccomplex_t IT_0209 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0210 = IT_0020*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*(IT_0204 + IT_0206 + 
      -IT_0208 + -IT_0210);
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = m_N_3*IT_0212;
    const ccomplex_t IT_0214 = IT_0202 + IT_0213;
    const ccomplex_t IT_0215 = cpow((-2)*s_23 + IT_0079 + IT_0080 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0216 = IT_0153*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = IT_0214*IT_0217;
    const ccomplex_t IT_0219 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0220 = IT_0002*IT_0219;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = m_b*N_d1*e_em*IT_0061*U_sb_00;
    const ccomplex_t IT_0223 = IT_0060*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*(IT_0221 + 1.5*IT_0224);
    const ccomplex_t IT_0226 = (-0.333333333333333)*IT_0225;
    const ccomplex_t IT_0227 = cpow((-2)*s_23 + IT_0079 + IT_0080 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0228 = IT_0226*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0231 = IT_0003*IT_0230;
    const ccomplex_t IT_0232 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0233 = IT_0007*IT_0232;
    const ccomplex_t IT_0234 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0235 = IT_0003*IT_0234;
    const ccomplex_t IT_0236 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0237 = IT_0007*IT_0236;
    const ccomplex_t IT_0238 = IT_0231 + IT_0233 + IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0240 = IT_0017*IT_0239;
    const ccomplex_t IT_0241 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0242 = IT_0020*IT_0241;
    const ccomplex_t IT_0243 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0244 = IT_0017*IT_0243;
    const ccomplex_t IT_0245 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0246 = IT_0020*IT_0245;
    const ccomplex_t IT_0247 = -IT_0240 + -IT_0242 + -IT_0244 + -IT_0246;
    const ccomplex_t IT_0248 = IT_0238 + IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = (-0.5)*IT_0249;
    const ccomplex_t IT_0251 = m_N_1*IT_0250;
    const ccomplex_t IT_0252 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0253 = IT_0003*IT_0252;
    const ccomplex_t IT_0254 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0255 = IT_0007*IT_0254;
    const ccomplex_t IT_0256 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0257 = IT_0003*IT_0256;
    const ccomplex_t IT_0258 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0259 = IT_0007*IT_0258;
    const ccomplex_t IT_0260 = IT_0253 + IT_0255 + IT_0257 + IT_0259;
    const ccomplex_t IT_0261 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0262 = IT_0017*IT_0261;
    const ccomplex_t IT_0263 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0264 = IT_0020*IT_0263;
    const ccomplex_t IT_0265 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0266 = IT_0017*IT_0265;
    const ccomplex_t IT_0267 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0268 = IT_0020*IT_0267;
    const ccomplex_t IT_0269 = -IT_0262 + -IT_0264 + -IT_0266 + -IT_0268;
    const ccomplex_t IT_0270 = IT_0260 + IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = (-0.5)*IT_0271;
    const ccomplex_t IT_0273 = m_N_3*IT_0272;
    const ccomplex_t IT_0274 = IT_0251 + IT_0273;
    const ccomplex_t IT_0275 = IT_0229*IT_0274;
    const ccomplex_t IT_0276 = IT_0084 + IT_0141 + IT_0166 + IT_0184 + (-2)
      *IT_0191 + IT_0218 + IT_0275;
    const ccomplex_t IT_0277 = m_b*conj(N_d3)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0278 = IT_0060*IT_0277;
    const ccomplex_t IT_0279 = 1.4142135623731*IT_0278;
    const ccomplex_t IT_0280 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0281 = IT_0002*IT_0280;
    const ccomplex_t IT_0282 = 1.4142135623731*IT_0281;
    const ccomplex_t IT_0283 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0284 = IT_0016*IT_0283;
    const ccomplex_t IT_0285 = 1.4142135623731*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*(IT_0279 + 0.333333333333333
      *IT_0282 + -IT_0285);
    const ccomplex_t IT_0287 = 3*IT_0286;
    const ccomplex_t IT_0288 = 0.333333333333333*IT_0287;
    const ccomplex_t IT_0289 = 3*IT_0288;
    const ccomplex_t IT_0290 = (-0.333333333333333)*IT_0289;
    const ccomplex_t IT_0291 = (-3)*IT_0290;
    const ccomplex_t IT_0292 = 0.333333333333333*IT_0291;
    const ccomplex_t IT_0293 = 3*IT_0292;
    const ccomplex_t IT_0294 = 0.166666666666667*IT_0293;
    const ccomplex_t IT_0295 = IT_0187*IT_0294;
    const ccomplex_t IT_0296 = IT_0185*IT_0295;
    const ccomplex_t IT_0297 = (0 + _Complex_I*1)*IT_0296;
    const ccomplex_t IT_0298 = m_b*IT_0297;
    const ccomplex_t IT_0299 = m_b*conj(N_d2)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0300 = IT_0060*IT_0299;
    const ccomplex_t IT_0301 = 1.4142135623731*IT_0300;
    const ccomplex_t IT_0302 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0303 = IT_0002*IT_0302;
    const ccomplex_t IT_0304 = 1.4142135623731*IT_0303;
    const ccomplex_t IT_0305 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0306 = IT_0016*IT_0305;
    const ccomplex_t IT_0307 = 1.4142135623731*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*(IT_0301 + 0.333333333333333
      *IT_0304 + -IT_0307);
    const ccomplex_t IT_0309 = 3*IT_0308;
    const ccomplex_t IT_0310 = 0.333333333333333*IT_0309;
    const ccomplex_t IT_0311 = -IT_0310;
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = -IT_0312;
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = 3*IT_0314;
    const ccomplex_t IT_0316 = 0.333333333333333*IT_0315;
    const ccomplex_t IT_0317 = 3*IT_0316;
    const ccomplex_t IT_0318 = (-0.333333333333333)*IT_0317;
    const ccomplex_t IT_0319 = -IT_0318;
    const ccomplex_t IT_0320 = 3*IT_0319;
    const ccomplex_t IT_0321 = 0.166666666666667*IT_0320;
    const ccomplex_t IT_0322 = IT_0081*IT_0321;
    const ccomplex_t IT_0323 = (0 + _Complex_I*1)*IT_0322;
    const ccomplex_t IT_0324 = m_N_2*IT_0052;
    const ccomplex_t IT_0325 = m_N_3*IT_0030;
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = IT_0323*IT_0326;
    const ccomplex_t IT_0328 = 2*U_sb_00*((IT_0006*IT_0016*IT_0059*IT_0061
      *IT_0079 + (-0.5)*m_W*IT_0159*(IT_0016 + 1./3*IT_0015*IT_0160))*conj
      (U_sb_00) + 1./2*IT_0000*IT_0016*IT_0059*IT_0155*conj(U_sb_10)) + U_sb_10*
      (IT_0000*IT_0016*IT_0059*IT_0155*conj(U_sb_00) + 2*(IT_0006*IT_0016
      *IT_0059*IT_0061*IT_0079 + -1./3*m_W*IT_0015*IT_0159*IT_0160)*conj(U_sb_10
      ));
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*(e_em*IT_0328 +
       1.4142135623731*A_b*IT_0006*IT_0154*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0330 = 0.5*IT_0329;
    const ccomplex_t IT_0331 = IT_0294*IT_0330;
    const ccomplex_t IT_0332 = IT_0143*IT_0331;
    const ccomplex_t IT_0333 = (0 + _Complex_I*1)*IT_0332;
    const ccomplex_t IT_0334 = m_b*conj(N_d3)*e_em*IT_0061*U_sb_11;
    const ccomplex_t IT_0335 = IT_0060*IT_0334;
    const ccomplex_t IT_0336 = 1.4142135623731*IT_0335;
    const ccomplex_t IT_0337 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0338 = IT_0016*IT_0337;
    const ccomplex_t IT_0339 = 1.4142135623731*IT_0338;
    const ccomplex_t IT_0340 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0341 = IT_0002*IT_0340;
    const ccomplex_t IT_0342 = 1.4142135623731*IT_0341;
    const ccomplex_t IT_0343 = (0 + _Complex_I*1)*(IT_0336 + -IT_0339 +
       0.333333333333333*IT_0342);
    const ccomplex_t IT_0344 = 3*IT_0343;
    const ccomplex_t IT_0345 = 0.166666666666667*IT_0344;
    const ccomplex_t IT_0346 = 2*U_sb_00*((IT_0006*IT_0016*IT_0059*IT_0061
      *IT_0079 + (-0.5)*m_W*IT_0159*(IT_0016 + 1./3*IT_0015*IT_0160))*conj
      (U_sb_01) + 1./2*IT_0000*IT_0016*IT_0059*IT_0155*conj(U_sb_11)) + U_sb_10*
      (IT_0000*IT_0016*IT_0059*IT_0155*conj(U_sb_01) + 2*(IT_0006*IT_0016
      *IT_0059*IT_0061*IT_0079 + -1./3*m_W*IT_0015*IT_0159*IT_0160)*conj(U_sb_11
      ));
    const ccomplex_t IT_0347 = (0 + _Complex_I*1)*(e_em*IT_0346 +
       1.4142135623731*A_b*IT_0006*IT_0154*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0348 = 0.5*IT_0347;
    const ccomplex_t IT_0349 = IT_0345*IT_0348;
    const ccomplex_t IT_0350 = IT_0167*IT_0349;
    const ccomplex_t IT_0351 = (0 + _Complex_I*1)*IT_0350;
    const ccomplex_t IT_0352 = IT_0215*IT_0294;
    const ccomplex_t IT_0353 = (0 + _Complex_I*1)*IT_0352;
    const ccomplex_t IT_0354 = IT_0214*IT_0353;
    const ccomplex_t IT_0355 = m_N_3*IT_0250;
    const ccomplex_t IT_0356 = m_N_1*IT_0272;
    const ccomplex_t IT_0357 = IT_0355 + IT_0356;
    const ccomplex_t IT_0358 = m_b*conj(N_d1)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0359 = IT_0060*IT_0358;
    const ccomplex_t IT_0360 = 1.4142135623731*IT_0359;
    const ccomplex_t IT_0361 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0362 = IT_0002*IT_0361;
    const ccomplex_t IT_0363 = 1.4142135623731*IT_0362;
    const ccomplex_t IT_0364 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0365 = IT_0016*IT_0364;
    const ccomplex_t IT_0366 = 1.4142135623731*IT_0365;
    const ccomplex_t IT_0367 = (0 + _Complex_I*1)*(IT_0360 + 0.333333333333333
      *IT_0363 + -IT_0366);
    const ccomplex_t IT_0368 = 3*IT_0367;
    const ccomplex_t IT_0369 = (-0.333333333333333)*IT_0368;
    const ccomplex_t IT_0370 = -IT_0369;
    const ccomplex_t IT_0371 = -IT_0370;
    const ccomplex_t IT_0372 = (-3)*IT_0371;
    const ccomplex_t IT_0373 = (-0.333333333333333)*IT_0372;
    const ccomplex_t IT_0374 = -IT_0373;
    const ccomplex_t IT_0375 = 3*IT_0374;
    const ccomplex_t IT_0376 = 0.166666666666667*IT_0375;
    const ccomplex_t IT_0377 = IT_0227*IT_0376;
    const ccomplex_t IT_0378 = (0 + _Complex_I*1)*IT_0377;
    const ccomplex_t IT_0379 = IT_0357*IT_0378;
    const ccomplex_t IT_0380 = m_b*conj(N_d4)*e_em*IT_0061*U_sb_10;
    const ccomplex_t IT_0381 = IT_0060*IT_0380;
    const ccomplex_t IT_0382 = 1.4142135623731*IT_0381;
    const ccomplex_t IT_0383 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0384 = IT_0016*IT_0383;
    const ccomplex_t IT_0385 = 1.4142135623731*IT_0384;
    const ccomplex_t IT_0386 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0387 = IT_0002*IT_0386;
    const ccomplex_t IT_0388 = 1.4142135623731*IT_0387;
    const ccomplex_t IT_0389 = (0 + _Complex_I*1)*(IT_0382 + -IT_0385 +
       0.333333333333333*IT_0388);
    const ccomplex_t IT_0390 = 3*IT_0389;
    const ccomplex_t IT_0391 = 0.166666666666667*IT_0390;
    const ccomplex_t IT_0392 = IT_0138*IT_0391;
    const ccomplex_t IT_0393 = (0 + _Complex_I*1)*IT_0392;
    const ccomplex_t IT_0394 = m_N_4*IT_0127;
    const ccomplex_t IT_0395 = m_N_3*IT_0105;
    const ccomplex_t IT_0396 = IT_0394 + IT_0395;
    const ccomplex_t IT_0397 = IT_0393*IT_0396;
    const ccomplex_t IT_0398 = 2*IT_0298 + -IT_0327 + -IT_0333 + -IT_0351 + 
      -IT_0354 + -IT_0379 + -IT_0397;
    const ccomplex_t IT_0399 = 6*s_13;
    const ccomplex_t IT_0400 = m_b*m_N_3;
    const ccomplex_t IT_0401 = 6*IT_0276;
    const ccomplex_t IT_0402 = IT_0250*IT_0378;
    const ccomplex_t IT_0403 = IT_0105*IT_0393;
    const ccomplex_t IT_0404 = IT_0030*IT_0323;
    const ccomplex_t IT_0405 = IT_0201*IT_0353;
    const ccomplex_t IT_0406 = -IT_0190 + IT_0402 + IT_0403 + IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = m_b*s_14;
    const ccomplex_t IT_0408 = IT_0127*IT_0140;
    const ccomplex_t IT_0409 = IT_0052*IT_0083;
    const ccomplex_t IT_0410 = IT_0212*IT_0217;
    const ccomplex_t IT_0411 = IT_0229*IT_0272;
    const ccomplex_t IT_0412 = IT_0297 + -IT_0408 + -IT_0409 + -IT_0410 + 
      -IT_0411;
    const ccomplex_t IT_0413 = s_34*m_N_3;
    const ccomplex_t IT_0414 = 6*conj(IT_0276);
    const ccomplex_t IT_0415 = 6*conj(IT_0398);
    const ccomplex_t IT_0416 = pow(m_h, 2);
    const ccomplex_t IT_0417 = m_b*m_N_3*IT_0416;
    const ccomplex_t IT_0418 = s_13*IT_0416;
    const ccomplex_t IT_0419 = s_14*s_34;
    const ccomplex_t IT_0420 = (-2)*IT_0419;
    const ccomplex_t IT_0421 = IT_0418 + IT_0420;
    const ccomplex_t IT_0422 = 6*IT_0398;
    const ccomplex_t IT_0423 = 0.166666666666667*IT_0413;
    const ccomplex_t IT_0424 = (IT_0276*conj(IT_0276) + IT_0398*conj(IT_0398))
      *IT_0399 + conj(IT_0398)*IT_0400*IT_0401 + (IT_0398*IT_0400 + IT_0406
      *IT_0407 + IT_0412*IT_0413)*IT_0414 + (IT_0407*IT_0412 + IT_0406*IT_0413)
      *IT_0415 + 6*conj(IT_0412)*(IT_0406*IT_0417 + -IT_0412*IT_0421 +
       0.166666666666667*IT_0407*IT_0422 + IT_0401*IT_0423) + 6*conj(IT_0406)*
      (0.166666666666667*IT_0401*IT_0407 + IT_0412*IT_0417 + -IT_0406*IT_0421 +
       IT_0422*IT_0423);
    return create_ccomplex_return(IT_0424);
}

