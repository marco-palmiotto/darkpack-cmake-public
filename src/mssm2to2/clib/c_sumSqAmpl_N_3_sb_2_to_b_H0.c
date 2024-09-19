#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_sb_2_to_b_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_sb_2_to_b_H0(
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = m_b*N_d2*e_em*IT_0010*U_sb_01;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0004 + 1.5*IT_0013);
    const ccomplex_t IT_0015 = 0.666666666666667*IT_0014;
    const ccomplex_t IT_0016 = 1.5*IT_0015;
    const ccomplex_t IT_0017 = (-0.333333333333333)*IT_0016;
    const ccomplex_t IT_0018 = pow(m_b, 2);
    const ccomplex_t IT_0019 = pow(m_sb_2, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = cos(alpha);
    const ccomplex_t IT_0024 = IT_0001*IT_0023;
    const ccomplex_t IT_0025 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = sin(alpha);
    const ccomplex_t IT_0028 = IT_0008*IT_0027;
    const ccomplex_t IT_0029 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0032 = IT_0024*IT_0031;
    const ccomplex_t IT_0033 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = IT_0026 + IT_0030 + IT_0032 + IT_0034;
    const ccomplex_t IT_0036 = IT_0001*IT_0027;
    const ccomplex_t IT_0037 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = IT_0008*IT_0023;
    const ccomplex_t IT_0040 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0043 = IT_0036*IT_0042;
    const ccomplex_t IT_0044 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0045 = IT_0039*IT_0044;
    const ccomplex_t IT_0046 = -IT_0038 + -IT_0041 + -IT_0043 + -IT_0045;
    const ccomplex_t IT_0047 = IT_0035 + IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = m_N_2*IT_0049;
    const ccomplex_t IT_0051 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0052 = IT_0024*IT_0051;
    const ccomplex_t IT_0053 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0054 = IT_0028*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0056 = IT_0024*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0058 = IT_0028*IT_0057;
    const ccomplex_t IT_0059 = IT_0052 + IT_0054 + IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0061 = IT_0036*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0063 = IT_0039*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0065 = IT_0036*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0067 = IT_0039*IT_0066;
    const ccomplex_t IT_0068 = -IT_0061 + -IT_0063 + -IT_0065 + -IT_0067;
    const ccomplex_t IT_0069 = IT_0059 + IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = m_N_3*IT_0071;
    const ccomplex_t IT_0073 = IT_0050 + IT_0072;
    const ccomplex_t IT_0074 = IT_0022*IT_0073;
    const ccomplex_t IT_0075 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0076 = IT_0001*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = m_b*N_d4*e_em*IT_0010*U_sb_01;
    const ccomplex_t IT_0079 = IT_0009*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 1.5*IT_0080);
    const ccomplex_t IT_0082 = 0.666666666666667*IT_0081;
    const ccomplex_t IT_0083 = 1.5*IT_0082;
    const ccomplex_t IT_0084 = (-0.333333333333333)*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0089 = IT_0024*IT_0088;
    const ccomplex_t IT_0090 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0091 = IT_0028*IT_0090;
    const ccomplex_t IT_0092 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0093 = IT_0024*IT_0092;
    const ccomplex_t IT_0094 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0095 = IT_0028*IT_0094;
    const ccomplex_t IT_0096 = IT_0089 + IT_0091 + IT_0093 + IT_0095;
    const ccomplex_t IT_0097 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0098 = IT_0036*IT_0097;
    const ccomplex_t IT_0099 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0100 = IT_0039*IT_0099;
    const ccomplex_t IT_0101 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0102 = IT_0036*IT_0101;
    const ccomplex_t IT_0103 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0104 = IT_0039*IT_0103;
    const ccomplex_t IT_0105 = -IT_0098 + -IT_0100 + -IT_0102 + -IT_0104;
    const ccomplex_t IT_0106 = IT_0096 + IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = m_N_4*IT_0108;
    const ccomplex_t IT_0110 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0111 = IT_0024*IT_0110;
    const ccomplex_t IT_0112 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0113 = IT_0028*IT_0112;
    const ccomplex_t IT_0114 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0115 = IT_0024*IT_0114;
    const ccomplex_t IT_0116 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0117 = IT_0028*IT_0116;
    const ccomplex_t IT_0118 = IT_0111 + IT_0113 + IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0120 = IT_0036*IT_0119;
    const ccomplex_t IT_0121 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0122 = IT_0039*IT_0121;
    const ccomplex_t IT_0123 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0124 = IT_0036*IT_0123;
    const ccomplex_t IT_0125 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0126 = IT_0039*IT_0125;
    const ccomplex_t IT_0127 = -IT_0120 + -IT_0122 + -IT_0124 + -IT_0126;
    const ccomplex_t IT_0128 = IT_0118 + IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = m_N_3*IT_0130;
    const ccomplex_t IT_0132 = IT_0109 + IT_0131;
    const ccomplex_t IT_0133 = IT_0087*IT_0132;
    const ccomplex_t IT_0134 = pow(m_N_3, 2);
    const ccomplex_t IT_0135 = cpow((-2)*s_13 + IT_0018 + IT_0134 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0136 = cpow(V_tb, 2);
    const ccomplex_t IT_0137 = m_b*mu_h*IT_0010;
    const ccomplex_t IT_0138 = sin(beta);
    const ccomplex_t IT_0139 = -IT_0027*IT_0138;
    const ccomplex_t IT_0140 = IT_0005*IT_0023;
    const ccomplex_t IT_0141 = IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = cpow(IT_0000, -2);
    const ccomplex_t IT_0143 = (-2)*U_sb_01*((IT_0006*IT_0008*IT_0010*IT_0018
      *IT_0023 + (-0.5)*m_W*IT_0141*(IT_0008 + 1./3*IT_0007*IT_0142))*conj
      (U_sb_01) + -1./2*IT_0006*IT_0008*IT_0027*IT_0137*conj(U_sb_11)) + U_sb_11
      *(IT_0006*IT_0008*IT_0027*IT_0137*conj(U_sb_01) + (-2)*(IT_0006*IT_0008
      *IT_0010*IT_0018*IT_0023 + -1./3*m_W*IT_0007*IT_0141*IT_0142)*conj(U_sb_11
      ));
    const ccomplex_t IT_0144 = (0 + _Complex_I*-1)*(e_em*IT_0143 + (
      -1.4142135623731)*A_b*IT_0023*IT_0136*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0145 = (-0.5)*IT_0144;
    const ccomplex_t IT_0146 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0147 = IT_0001*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = m_b*N_d3*e_em*IT_0010*U_sb_01;
    const ccomplex_t IT_0150 = IT_0009*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0148 + 1.5*IT_0151);
    const ccomplex_t IT_0153 = 0.666666666666667*IT_0152;
    const ccomplex_t IT_0154 = 1.5*IT_0153;
    const ccomplex_t IT_0155 = 0.666666666666667*IT_0154;
    const ccomplex_t IT_0156 = 1.5*IT_0155;
    const ccomplex_t IT_0157 = 0.666666666666667*IT_0156;
    const ccomplex_t IT_0158 = 1.5*IT_0157;
    const ccomplex_t IT_0159 = (-0.333333333333333)*IT_0158;
    const ccomplex_t IT_0160 = IT_0145*IT_0159;
    const ccomplex_t IT_0161 = IT_0135*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = cpow((-2)*s_13 + IT_0018 + IT_0134 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0164 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0165 = IT_0001*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = m_b*N_d3*e_em*IT_0010*U_sb_00;
    const ccomplex_t IT_0168 = IT_0009*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*(IT_0166 + 1.5*IT_0169);
    const ccomplex_t IT_0171 = 0.666666666666667*IT_0170;
    const ccomplex_t IT_0172 = 1.5*IT_0171;
    const ccomplex_t IT_0173 = (-0.333333333333333)*IT_0172;
    const ccomplex_t IT_0174 = (-2)*conj(U_sb_00)*((IT_0006*IT_0008*IT_0010
      *IT_0018*IT_0023 + (-0.5)*m_W*IT_0141*(IT_0008 + 1./3*IT_0007*IT_0142))
      *U_sb_01 + -1./2*IT_0006*IT_0008*IT_0027*IT_0137*U_sb_11) + conj(U_sb_10)*
      (IT_0006*IT_0008*IT_0027*IT_0137*U_sb_01 + (-2)*(IT_0006*IT_0008*IT_0010
      *IT_0018*IT_0023 + -1./3*m_W*IT_0007*IT_0141*IT_0142)*U_sb_11);
    const ccomplex_t IT_0175 = (0 + _Complex_I*-1)*(e_em*IT_0174 + (
      -1.4142135623731)*A_b*IT_0023*IT_0136*(U_sb_01*conj(U_sb_10) + conj
      (U_sb_00)*U_sb_11));
    const ccomplex_t IT_0176 = (-0.5)*IT_0175;
    const ccomplex_t IT_0177 = IT_0173*IT_0176;
    const ccomplex_t IT_0178 = IT_0163*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = cpow((-2)*s_12 + IT_0018 + -IT_0019 + -IT_0134 
      + -reg_prop, -1);
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*m_b*e_em*IT_0006*IT_0008
      *IT_0010*IT_0023;
    const ccomplex_t IT_0182 = (-0.5)*IT_0181;
    const ccomplex_t IT_0183 = IT_0159*IT_0182;
    const ccomplex_t IT_0184 = IT_0180*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = m_b*IT_0185;
    const ccomplex_t IT_0187 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0188 = IT_0024*IT_0187;
    const ccomplex_t IT_0189 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0190 = IT_0036*IT_0189;
    const ccomplex_t IT_0191 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0192 = IT_0039*IT_0191;
    const ccomplex_t IT_0193 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0194 = IT_0028*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*(IT_0188 + -IT_0190 + 
      -IT_0192 + IT_0194);
    const ccomplex_t IT_0196 = m_N_3*IT_0195;
    const ccomplex_t IT_0197 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0198 = IT_0024*IT_0197;
    const ccomplex_t IT_0199 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0200 = IT_0036*IT_0199;
    const ccomplex_t IT_0201 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0202 = IT_0039*IT_0201;
    const ccomplex_t IT_0203 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0204 = IT_0028*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*(IT_0198 + -IT_0200 + 
      -IT_0202 + IT_0204);
    const ccomplex_t IT_0206 = m_N_3*IT_0205;
    const ccomplex_t IT_0207 = IT_0196 + IT_0206;
    const ccomplex_t IT_0208 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0209 = IT_0159*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = IT_0207*IT_0210;
    const ccomplex_t IT_0212 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0213 = IT_0024*IT_0212;
    const ccomplex_t IT_0214 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0215 = IT_0028*IT_0214;
    const ccomplex_t IT_0216 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0217 = IT_0024*IT_0216;
    const ccomplex_t IT_0218 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0219 = IT_0028*IT_0218;
    const ccomplex_t IT_0220 = IT_0213 + IT_0215 + IT_0217 + IT_0219;
    const ccomplex_t IT_0221 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0222 = IT_0036*IT_0221;
    const ccomplex_t IT_0223 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0224 = IT_0039*IT_0223;
    const ccomplex_t IT_0225 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0226 = IT_0036*IT_0225;
    const ccomplex_t IT_0227 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0228 = IT_0039*IT_0227;
    const ccomplex_t IT_0229 = -IT_0222 + -IT_0224 + -IT_0226 + -IT_0228;
    const ccomplex_t IT_0230 = IT_0220 + IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = 0.5*IT_0231;
    const ccomplex_t IT_0233 = m_N_1*IT_0232;
    const ccomplex_t IT_0234 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0235 = IT_0024*IT_0234;
    const ccomplex_t IT_0236 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0237 = IT_0028*IT_0236;
    const ccomplex_t IT_0238 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0239 = IT_0024*IT_0238;
    const ccomplex_t IT_0240 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0241 = IT_0028*IT_0240;
    const ccomplex_t IT_0242 = IT_0235 + IT_0237 + IT_0239 + IT_0241;
    const ccomplex_t IT_0243 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0244 = IT_0036*IT_0243;
    const ccomplex_t IT_0245 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0246 = IT_0039*IT_0245;
    const ccomplex_t IT_0247 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0248 = IT_0036*IT_0247;
    const ccomplex_t IT_0249 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0250 = IT_0039*IT_0249;
    const ccomplex_t IT_0251 = -IT_0244 + -IT_0246 + -IT_0248 + -IT_0250;
    const ccomplex_t IT_0252 = IT_0242 + IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = 0.5*IT_0253;
    const ccomplex_t IT_0255 = m_N_3*IT_0254;
    const ccomplex_t IT_0256 = IT_0233 + IT_0255;
    const ccomplex_t IT_0257 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0258 = IT_0001*IT_0257;
    const ccomplex_t IT_0259 = 1.4142135623731*IT_0258;
    const ccomplex_t IT_0260 = m_b*N_d1*e_em*IT_0010*U_sb_01;
    const ccomplex_t IT_0261 = IT_0009*IT_0260;
    const ccomplex_t IT_0262 = 1.4142135623731*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*(IT_0259 + 1.5*IT_0262);
    const ccomplex_t IT_0264 = (-0.333333333333333)*IT_0263;
    const ccomplex_t IT_0265 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0266 = IT_0264*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = IT_0256*IT_0267;
    const ccomplex_t IT_0269 = IT_0074 + IT_0133 + IT_0162 + IT_0179 + (-2)
      *IT_0186 + IT_0211 + IT_0268;
    const ccomplex_t IT_0270 = m_b*conj(N_d2)*e_em*IT_0010*U_sb_11;
    const ccomplex_t IT_0271 = IT_0009*IT_0270;
    const ccomplex_t IT_0272 = 1.4142135623731*IT_0271;
    const ccomplex_t IT_0273 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0274 = IT_0001*IT_0273;
    const ccomplex_t IT_0275 = 1.4142135623731*IT_0274;
    const ccomplex_t IT_0276 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0277 = IT_0008*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*(IT_0272 + 0.333333333333333
      *IT_0275 + -IT_0278);
    const ccomplex_t IT_0280 = -IT_0279;
    const ccomplex_t IT_0281 = (-3)*IT_0280;
    const ccomplex_t IT_0282 = 0.166666666666667*IT_0281;
    const ccomplex_t IT_0283 = IT_0020*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = m_N_3*IT_0049;
    const ccomplex_t IT_0286 = m_N_2*IT_0071;
    const ccomplex_t IT_0287 = IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = IT_0284*IT_0287;
    const ccomplex_t IT_0289 = m_b*conj(N_d3)*e_em*IT_0010*U_sb_11;
    const ccomplex_t IT_0290 = IT_0009*IT_0289;
    const ccomplex_t IT_0291 = 1.4142135623731*IT_0290;
    const ccomplex_t IT_0292 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0293 = IT_0008*IT_0292;
    const ccomplex_t IT_0294 = 1.4142135623731*IT_0293;
    const ccomplex_t IT_0295 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0296 = IT_0001*IT_0295;
    const ccomplex_t IT_0297 = 1.4142135623731*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*(IT_0291 + -IT_0294 +
       0.333333333333333*IT_0297);
    const ccomplex_t IT_0299 = 3*IT_0298;
    const ccomplex_t IT_0300 = 0.166666666666667*IT_0299;
    const ccomplex_t IT_0301 = IT_0208*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*IT_0301;
    const ccomplex_t IT_0303 = IT_0207*IT_0302;
    const ccomplex_t IT_0304 = m_b*conj(N_d1)*e_em*IT_0010*U_sb_11;
    const ccomplex_t IT_0305 = IT_0009*IT_0304;
    const ccomplex_t IT_0306 = 1.4142135623731*IT_0305;
    const ccomplex_t IT_0307 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0308 = IT_0008*IT_0307;
    const ccomplex_t IT_0309 = 1.4142135623731*IT_0308;
    const ccomplex_t IT_0310 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0311 = IT_0001*IT_0310;
    const ccomplex_t IT_0312 = 1.4142135623731*IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*(IT_0306 + -IT_0309 +
       0.333333333333333*IT_0312);
    const ccomplex_t IT_0314 = 3*IT_0313;
    const ccomplex_t IT_0315 = 0.166666666666667*IT_0314;
    const ccomplex_t IT_0316 = IT_0265*IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*IT_0316;
    const ccomplex_t IT_0318 = m_N_3*IT_0232;
    const ccomplex_t IT_0319 = m_N_1*IT_0254;
    const ccomplex_t IT_0320 = IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = IT_0317*IT_0320;
    const ccomplex_t IT_0322 = m_b*conj(N_d4)*e_em*IT_0010*U_sb_11;
    const ccomplex_t IT_0323 = IT_0009*IT_0322;
    const ccomplex_t IT_0324 = 1.4142135623731*IT_0323;
    const ccomplex_t IT_0325 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0326 = IT_0001*IT_0325;
    const ccomplex_t IT_0327 = 1.4142135623731*IT_0326;
    const ccomplex_t IT_0328 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0329 = IT_0008*IT_0328;
    const ccomplex_t IT_0330 = 1.4142135623731*IT_0329;
    const ccomplex_t IT_0331 = (0 + _Complex_I*1)*(IT_0324 + 0.333333333333333
      *IT_0327 + -IT_0330);
    const ccomplex_t IT_0332 = 3*IT_0331;
    const ccomplex_t IT_0333 = 0.166666666666667*IT_0332;
    const ccomplex_t IT_0334 = IT_0085*IT_0333;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = m_N_3*IT_0108;
    const ccomplex_t IT_0337 = m_N_4*IT_0130;
    const ccomplex_t IT_0338 = IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = IT_0335*IT_0338;
    const ccomplex_t IT_0340 = (-2)*U_sb_01*((IT_0006*IT_0008*IT_0010*IT_0018
      *IT_0023 + (-0.5)*m_W*IT_0141*(IT_0008 + 1./3*IT_0007*IT_0142))*conj
      (U_sb_01) + -1./2*IT_0006*IT_0008*IT_0027*IT_0137*conj(U_sb_11)) + U_sb_11
      *(IT_0006*IT_0008*IT_0027*IT_0137*conj(U_sb_01) + (-2)*(IT_0006*IT_0008
      *IT_0010*IT_0018*IT_0023 + -1./3*m_W*IT_0007*IT_0141*IT_0142)*conj(U_sb_11
      ));
    const ccomplex_t IT_0341 = (0 + _Complex_I*-1)*(e_em*IT_0340 + (
      -1.4142135623731)*A_b*IT_0023*IT_0136*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0342 = (-0.5)*IT_0341;
    const ccomplex_t IT_0343 = IT_0300*IT_0342;
    const ccomplex_t IT_0344 = IT_0135*IT_0343;
    const ccomplex_t IT_0345 = (0 + _Complex_I*1)*IT_0344;
    const ccomplex_t IT_0346 = m_b*conj(N_d3)*e_em*IT_0010*U_sb_10;
    const ccomplex_t IT_0347 = IT_0009*IT_0346;
    const ccomplex_t IT_0348 = 1.4142135623731*IT_0347;
    const ccomplex_t IT_0349 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0350 = IT_0008*IT_0349;
    const ccomplex_t IT_0351 = 1.4142135623731*IT_0350;
    const ccomplex_t IT_0352 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0353 = IT_0001*IT_0352;
    const ccomplex_t IT_0354 = 1.4142135623731*IT_0353;
    const ccomplex_t IT_0355 = (0 + _Complex_I*1)*(IT_0348 + -IT_0351 +
       0.333333333333333*IT_0354);
    const ccomplex_t IT_0356 = 3*IT_0355;
    const ccomplex_t IT_0357 = (-0.333333333333333)*IT_0356;
    const ccomplex_t IT_0358 = -IT_0357;
    const ccomplex_t IT_0359 = 3*IT_0358;
    const ccomplex_t IT_0360 = (-0.333333333333333)*IT_0359;
    const ccomplex_t IT_0361 = -IT_0360;
    const ccomplex_t IT_0362 = -IT_0361;
    const ccomplex_t IT_0363 = -IT_0362;
    const ccomplex_t IT_0364 = -IT_0363;
    const ccomplex_t IT_0365 = -IT_0364;
    const ccomplex_t IT_0366 = 3*IT_0365;
    const ccomplex_t IT_0367 = (-0.333333333333333)*IT_0366;
    const ccomplex_t IT_0368 = (-3)*IT_0367;
    const ccomplex_t IT_0369 = 0.333333333333333*IT_0368;
    const ccomplex_t IT_0370 = 3*IT_0369;
    const ccomplex_t IT_0371 = 0.166666666666667*IT_0370;
    const ccomplex_t IT_0372 = (-2)*conj(U_sb_00)*((IT_0006*IT_0008*IT_0010
      *IT_0018*IT_0023 + (-0.5)*m_W*IT_0141*(IT_0008 + 1./3*IT_0007*IT_0142))
      *U_sb_01 + -1./2*IT_0006*IT_0008*IT_0027*IT_0137*U_sb_11) + conj(U_sb_10)*
      (IT_0006*IT_0008*IT_0027*IT_0137*U_sb_01 + (-2)*(IT_0006*IT_0008*IT_0010
      *IT_0018*IT_0023 + -1./3*m_W*IT_0007*IT_0141*IT_0142)*U_sb_11);
    const ccomplex_t IT_0373 = (0 + _Complex_I*-1)*(e_em*IT_0372 + (
      -1.4142135623731)*A_b*IT_0023*IT_0136*(U_sb_01*conj(U_sb_10) + conj
      (U_sb_00)*U_sb_11));
    const ccomplex_t IT_0374 = (-0.5)*IT_0373;
    const ccomplex_t IT_0375 = IT_0371*IT_0374;
    const ccomplex_t IT_0376 = IT_0163*IT_0375;
    const ccomplex_t IT_0377 = (0 + _Complex_I*1)*IT_0376;
    const ccomplex_t IT_0378 = IT_0182*IT_0300;
    const ccomplex_t IT_0379 = IT_0180*IT_0378;
    const ccomplex_t IT_0380 = (0 + _Complex_I*1)*IT_0379;
    const ccomplex_t IT_0381 = m_b*IT_0380;
    const ccomplex_t IT_0382 = -IT_0288 + -IT_0303 + -IT_0321 + -IT_0339 + 
      -IT_0345 + -IT_0377 + 2*IT_0381;
    const ccomplex_t IT_0383 = 6*s_13;
    const ccomplex_t IT_0384 = m_b*m_N_3;
    const ccomplex_t IT_0385 = 6*IT_0269;
    const ccomplex_t IT_0386 = IT_0195*IT_0302;
    const ccomplex_t IT_0387 = IT_0232*IT_0317;
    const ccomplex_t IT_0388 = IT_0049*IT_0284;
    const ccomplex_t IT_0389 = IT_0108*IT_0335;
    const ccomplex_t IT_0390 = -IT_0185 + IT_0386 + IT_0387 + IT_0388 + IT_0389;
    const ccomplex_t IT_0391 = m_b*s_14;
    const ccomplex_t IT_0392 = IT_0022*IT_0071;
    const ccomplex_t IT_0393 = IT_0254*IT_0267;
    const ccomplex_t IT_0394 = IT_0205*IT_0210;
    const ccomplex_t IT_0395 = IT_0087*IT_0130;
    const ccomplex_t IT_0396 = IT_0380 + -IT_0392 + -IT_0393 + -IT_0394 + 
      -IT_0395;
    const ccomplex_t IT_0397 = s_34*m_N_3;
    const ccomplex_t IT_0398 = 6*conj(IT_0269);
    const ccomplex_t IT_0399 = 6*conj(IT_0382);
    const ccomplex_t IT_0400 = pow(m_H0, 2);
    const ccomplex_t IT_0401 = m_b*m_N_3*IT_0400;
    const ccomplex_t IT_0402 = s_13*IT_0400;
    const ccomplex_t IT_0403 = s_14*s_34;
    const ccomplex_t IT_0404 = (-2)*IT_0403;
    const ccomplex_t IT_0405 = IT_0402 + IT_0404;
    const ccomplex_t IT_0406 = 6*IT_0382;
    const ccomplex_t IT_0407 = 0.166666666666667*IT_0397;
    const ccomplex_t IT_0408 = (IT_0269*conj(IT_0269) + IT_0382*conj(IT_0382))
      *IT_0383 + conj(IT_0382)*IT_0384*IT_0385 + (IT_0382*IT_0384 + IT_0390
      *IT_0391 + IT_0396*IT_0397)*IT_0398 + (IT_0391*IT_0396 + IT_0390*IT_0397)
      *IT_0399 + 6*conj(IT_0396)*(IT_0390*IT_0401 + -IT_0396*IT_0405 +
       0.166666666666667*IT_0391*IT_0406 + IT_0385*IT_0407) + 6*conj(IT_0390)*
      (0.166666666666667*IT_0385*IT_0391 + IT_0396*IT_0401 + -IT_0390*IT_0405 +
       IT_0406*IT_0407);
    return create_ccomplex_return(IT_0408);
}

