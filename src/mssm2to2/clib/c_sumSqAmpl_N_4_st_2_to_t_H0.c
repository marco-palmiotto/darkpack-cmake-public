#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_st_2_to_t_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_st_2_to_t_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
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
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
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
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = m_t*N_u4*e_em*IT_0010*U_st_01;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0004 + (-0.75)*IT_0013);
    const ccomplex_t IT_0015 = 0.666666666666667*IT_0014;
    const ccomplex_t IT_0016 = sin(alpha);
    const ccomplex_t IT_0017 = m_t*mu_h*IT_0010;
    const ccomplex_t IT_0018 = cos(alpha);
    const ccomplex_t IT_0019 = cos(beta);
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = -IT_0005*IT_0016 + IT_0020;
    const ccomplex_t IT_0022 = cpow(IT_0000, -2);
    const ccomplex_t IT_0023 = pow(m_t, 2);
    const ccomplex_t IT_0024 = (-2)*U_st_01*((0.5*m_W*IT_0021*(IT_0008 + -1./3
      *IT_0007*IT_0022) + IT_0006*IT_0008*IT_0010*IT_0016*IT_0023)*conj(U_st_01)
       + -1./2*IT_0006*IT_0008*IT_0017*IT_0018*conj(U_st_11)) + U_st_11*(IT_0006
      *IT_0008*IT_0017*IT_0018*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0007
      *IT_0021*IT_0022 + 1.5*IT_0006*IT_0008*IT_0010*IT_0016*IT_0023)*conj
      (U_st_11));
    const ccomplex_t IT_0025 = (0 + _Complex_I*-1)*(e_em*IT_0024 + (
      -1.4142135623731)*A_t*IT_0016*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0015*IT_0026;
    const ccomplex_t IT_0028 = pow(m_N_4, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0023 + IT_0028 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = m_t*N_u4*e_em*IT_0010*U_st_00;
    const ccomplex_t IT_0036 = IT_0009*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0034 + (-0.75)*IT_0037);
    const ccomplex_t IT_0039 = 0.666666666666667*IT_0038;
    const ccomplex_t IT_0040 = (-2)*conj(U_st_00)*((0.5*m_W*IT_0021*(IT_0008 +
       -1./3*IT_0007*IT_0022) + IT_0006*IT_0008*IT_0010*IT_0016*IT_0023)*U_st_01
       + -1./2*IT_0006*IT_0008*IT_0017*IT_0018*U_st_11) + conj(U_st_10)*(IT_0006
      *IT_0008*IT_0017*IT_0018*U_st_01 + (-1.33333333333333)*(m_W*IT_0007
      *IT_0021*IT_0022 + 1.5*IT_0006*IT_0008*IT_0010*IT_0016*IT_0023)*U_st_11);
    const ccomplex_t IT_0041 = (0 + _Complex_I*-1)*(e_em*IT_0040 + (
      -1.4142135623731)*A_t*IT_0016*(U_st_01*conj(U_st_10) + conj(U_st_00)
      *U_st_11));
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0039*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0023 + IT_0028 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*m_t*e_em*IT_0006*IT_0008
      *IT_0010*IT_0016;
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = IT_0015*IT_0048;
    const ccomplex_t IT_0050 = pow(m_st_2, 2);
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + -IT_0028 + -IT_0050 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = m_t*IT_0053;
    const ccomplex_t IT_0055 = IT_0001*IT_0018;
    const ccomplex_t IT_0056 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = IT_0001*IT_0016;
    const ccomplex_t IT_0059 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = IT_0008*IT_0018;
    const ccomplex_t IT_0062 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = IT_0008*IT_0016;
    const ccomplex_t IT_0065 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0057 + -IT_0060 + 
      -IT_0063 + IT_0066);
    const ccomplex_t IT_0068 = m_N_4*IT_0067;
    const ccomplex_t IT_0069 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0070 = IT_0055*IT_0069;
    const ccomplex_t IT_0071 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0072 = IT_0058*IT_0071;
    const ccomplex_t IT_0073 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0074 = IT_0061*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0076 = IT_0064*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0070 + -IT_0072 + 
      -IT_0074 + IT_0076);
    const ccomplex_t IT_0078 = m_N_4*IT_0077;
    const ccomplex_t IT_0079 = IT_0068 + IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_23 + IT_0023 + IT_0050 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0081 = IT_0015*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = IT_0079*IT_0082;
    const ccomplex_t IT_0084 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0085 = IT_0001*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = m_t*N_u1*e_em*IT_0010*U_st_01;
    const ccomplex_t IT_0088 = IT_0009*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0086 + (-0.75)*IT_0089);
    const ccomplex_t IT_0091 = 0.666666666666667*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_23 + IT_0023 + IT_0050 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0096 = IT_0055*IT_0095;
    const ccomplex_t IT_0097 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0098 = IT_0064*IT_0097;
    const ccomplex_t IT_0099 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0100 = IT_0055*IT_0099;
    const ccomplex_t IT_0101 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0102 = IT_0064*IT_0101;
    const ccomplex_t IT_0103 = IT_0096 + IT_0098 + IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0105 = IT_0058*IT_0104;
    const ccomplex_t IT_0106 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0107 = IT_0061*IT_0106;
    const ccomplex_t IT_0108 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0109 = IT_0058*IT_0108;
    const ccomplex_t IT_0110 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0111 = IT_0061*IT_0110;
    const ccomplex_t IT_0112 = -IT_0105 + -IT_0107 + -IT_0109 + -IT_0111;
    const ccomplex_t IT_0113 = IT_0103 + IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = m_N_1*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0118 = IT_0055*IT_0117;
    const ccomplex_t IT_0119 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0120 = IT_0064*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0122 = IT_0055*IT_0121;
    const ccomplex_t IT_0123 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0124 = IT_0064*IT_0123;
    const ccomplex_t IT_0125 = IT_0118 + IT_0120 + IT_0122 + IT_0124;
    const ccomplex_t IT_0126 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0127 = IT_0058*IT_0126;
    const ccomplex_t IT_0128 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0129 = IT_0061*IT_0128;
    const ccomplex_t IT_0130 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0131 = IT_0058*IT_0130;
    const ccomplex_t IT_0132 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0133 = IT_0061*IT_0132;
    const ccomplex_t IT_0134 = -IT_0127 + -IT_0129 + -IT_0131 + -IT_0133;
    const ccomplex_t IT_0135 = IT_0125 + IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = m_N_4*IT_0137;
    const ccomplex_t IT_0139 = IT_0116 + IT_0138;
    const ccomplex_t IT_0140 = IT_0094*IT_0139;
    const ccomplex_t IT_0141 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0142 = IT_0001*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = m_t*N_u2*e_em*IT_0010*U_st_01;
    const ccomplex_t IT_0145 = IT_0009*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0143 + (-0.75)*IT_0146);
    const ccomplex_t IT_0148 = 0.666666666666667*IT_0147;
    const ccomplex_t IT_0149 = cpow((-2)*s_23 + IT_0023 + IT_0050 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0150 = IT_0148*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0153 = IT_0055*IT_0152;
    const ccomplex_t IT_0154 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0155 = IT_0064*IT_0154;
    const ccomplex_t IT_0156 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0157 = IT_0055*IT_0156;
    const ccomplex_t IT_0158 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0159 = IT_0064*IT_0158;
    const ccomplex_t IT_0160 = IT_0153 + IT_0155 + IT_0157 + IT_0159;
    const ccomplex_t IT_0161 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0162 = IT_0058*IT_0161;
    const ccomplex_t IT_0163 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0164 = IT_0061*IT_0163;
    const ccomplex_t IT_0165 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0166 = IT_0058*IT_0165;
    const ccomplex_t IT_0167 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0168 = IT_0061*IT_0167;
    const ccomplex_t IT_0169 = -IT_0162 + -IT_0164 + -IT_0166 + -IT_0168;
    const ccomplex_t IT_0170 = IT_0160 + IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = 0.5*IT_0171;
    const ccomplex_t IT_0173 = m_N_4*IT_0172;
    const ccomplex_t IT_0174 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0175 = IT_0055*IT_0174;
    const ccomplex_t IT_0176 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0177 = IT_0064*IT_0176;
    const ccomplex_t IT_0178 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0179 = IT_0055*IT_0178;
    const ccomplex_t IT_0180 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0181 = IT_0064*IT_0180;
    const ccomplex_t IT_0182 = IT_0175 + IT_0177 + IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0184 = IT_0058*IT_0183;
    const ccomplex_t IT_0185 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0186 = IT_0061*IT_0185;
    const ccomplex_t IT_0187 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0188 = IT_0058*IT_0187;
    const ccomplex_t IT_0189 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0190 = IT_0061*IT_0189;
    const ccomplex_t IT_0191 = -IT_0184 + -IT_0186 + -IT_0188 + -IT_0190;
    const ccomplex_t IT_0192 = IT_0182 + IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = 0.5*IT_0193;
    const ccomplex_t IT_0195 = m_N_2*IT_0194;
    const ccomplex_t IT_0196 = IT_0173 + IT_0195;
    const ccomplex_t IT_0197 = IT_0151*IT_0196;
    const ccomplex_t IT_0198 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0199 = IT_0001*IT_0198;
    const ccomplex_t IT_0200 = 1.4142135623731*IT_0199;
    const ccomplex_t IT_0201 = m_t*N_u3*e_em*IT_0010*U_st_01;
    const ccomplex_t IT_0202 = IT_0009*IT_0201;
    const ccomplex_t IT_0203 = 1.4142135623731*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*(IT_0200 + (-0.75)*IT_0203);
    const ccomplex_t IT_0205 = 0.666666666666667*IT_0204;
    const ccomplex_t IT_0206 = cpow((-2)*s_23 + IT_0023 + IT_0050 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0207 = IT_0205*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0210 = IT_0055*IT_0209;
    const ccomplex_t IT_0211 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0212 = IT_0064*IT_0211;
    const ccomplex_t IT_0213 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0214 = IT_0055*IT_0213;
    const ccomplex_t IT_0215 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0216 = IT_0064*IT_0215;
    const ccomplex_t IT_0217 = IT_0210 + IT_0212 + IT_0214 + IT_0216;
    const ccomplex_t IT_0218 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0219 = IT_0058*IT_0218;
    const ccomplex_t IT_0220 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0221 = IT_0061*IT_0220;
    const ccomplex_t IT_0222 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0223 = IT_0058*IT_0222;
    const ccomplex_t IT_0224 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0225 = IT_0061*IT_0224;
    const ccomplex_t IT_0226 = -IT_0219 + -IT_0221 + -IT_0223 + -IT_0225;
    const ccomplex_t IT_0227 = IT_0217 + IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = 0.5*IT_0228;
    const ccomplex_t IT_0230 = m_N_3*IT_0229;
    const ccomplex_t IT_0231 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0232 = IT_0055*IT_0231;
    const ccomplex_t IT_0233 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0234 = IT_0064*IT_0233;
    const ccomplex_t IT_0235 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0236 = IT_0055*IT_0235;
    const ccomplex_t IT_0237 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0238 = IT_0064*IT_0237;
    const ccomplex_t IT_0239 = IT_0232 + IT_0234 + IT_0236 + IT_0238;
    const ccomplex_t IT_0240 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0241 = IT_0058*IT_0240;
    const ccomplex_t IT_0242 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0243 = IT_0061*IT_0242;
    const ccomplex_t IT_0244 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0245 = IT_0058*IT_0244;
    const ccomplex_t IT_0246 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0247 = IT_0061*IT_0246;
    const ccomplex_t IT_0248 = -IT_0241 + -IT_0243 + -IT_0245 + -IT_0247;
    const ccomplex_t IT_0249 = IT_0239 + IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*IT_0249;
    const ccomplex_t IT_0251 = 0.5*IT_0250;
    const ccomplex_t IT_0252 = m_N_4*IT_0251;
    const ccomplex_t IT_0253 = IT_0230 + IT_0252;
    const ccomplex_t IT_0254 = IT_0208*IT_0253;
    const ccomplex_t IT_0255 = IT_0031 + IT_0046 + (-2)*IT_0054 + IT_0083 +
       IT_0140 + IT_0197 + IT_0254;
    const ccomplex_t IT_0256 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0257 = IT_0001*IT_0256;
    const ccomplex_t IT_0258 = 1.4142135623731*IT_0257;
    const ccomplex_t IT_0259 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0260 = IT_0008*IT_0259;
    const ccomplex_t IT_0261 = 1.4142135623731*IT_0260;
    const ccomplex_t IT_0262 = m_t*conj(N_u4)*e_em*IT_0010*U_st_11;
    const ccomplex_t IT_0263 = IT_0009*IT_0262;
    const ccomplex_t IT_0264 = 1.4142135623731*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*(IT_0258 + 3*IT_0261 + 3
      *IT_0264);
    const ccomplex_t IT_0266 = 0.166666666666667*IT_0265;
    const ccomplex_t IT_0267 = IT_0026*IT_0266;
    const ccomplex_t IT_0268 = IT_0029*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*IT_0268;
    const ccomplex_t IT_0270 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0271 = IT_0001*IT_0270;
    const ccomplex_t IT_0272 = 1.4142135623731*IT_0271;
    const ccomplex_t IT_0273 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0274 = IT_0008*IT_0273;
    const ccomplex_t IT_0275 = 1.4142135623731*IT_0274;
    const ccomplex_t IT_0276 = m_t*conj(N_u4)*e_em*IT_0010*U_st_10;
    const ccomplex_t IT_0277 = IT_0009*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*(IT_0272 + 3*IT_0275 + 3
      *IT_0278);
    const ccomplex_t IT_0280 = 0.166666666666667*IT_0279;
    const ccomplex_t IT_0281 = IT_0042*IT_0280;
    const ccomplex_t IT_0282 = IT_0044*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*IT_0282;
    const ccomplex_t IT_0284 = IT_0048*IT_0266;
    const ccomplex_t IT_0285 = IT_0051*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = m_t*IT_0286;
    const ccomplex_t IT_0288 = IT_0080*IT_0266;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = IT_0079*IT_0289;
    const ccomplex_t IT_0291 = m_N_4*IT_0115;
    const ccomplex_t IT_0292 = m_N_1*IT_0137;
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0295 = IT_0001*IT_0294;
    const ccomplex_t IT_0296 = 1.4142135623731*IT_0295;
    const ccomplex_t IT_0297 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0298 = IT_0008*IT_0297;
    const ccomplex_t IT_0299 = 1.4142135623731*IT_0298;
    const ccomplex_t IT_0300 = m_t*conj(N_u1)*e_em*IT_0010*U_st_11;
    const ccomplex_t IT_0301 = IT_0009*IT_0300;
    const ccomplex_t IT_0302 = 1.4142135623731*IT_0301;
    const ccomplex_t IT_0303 = (0 + _Complex_I*1)*(IT_0296 + 3*IT_0299 + 3
      *IT_0302);
    const ccomplex_t IT_0304 = 0.166666666666667*IT_0303;
    const ccomplex_t IT_0305 = IT_0092*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*IT_0305;
    const ccomplex_t IT_0307 = IT_0293*IT_0306;
    const ccomplex_t IT_0308 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0309 = IT_0001*IT_0308;
    const ccomplex_t IT_0310 = 1.4142135623731*IT_0309;
    const ccomplex_t IT_0311 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0312 = IT_0008*IT_0311;
    const ccomplex_t IT_0313 = 1.4142135623731*IT_0312;
    const ccomplex_t IT_0314 = m_t*conj(N_u2)*e_em*IT_0010*U_st_11;
    const ccomplex_t IT_0315 = IT_0009*IT_0314;
    const ccomplex_t IT_0316 = 1.4142135623731*IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*(IT_0310 + 3*IT_0313 + 3
      *IT_0316);
    const ccomplex_t IT_0318 = 0.166666666666667*IT_0317;
    const ccomplex_t IT_0319 = IT_0149*IT_0318;
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*IT_0319;
    const ccomplex_t IT_0321 = m_N_2*IT_0172;
    const ccomplex_t IT_0322 = m_N_4*IT_0194;
    const ccomplex_t IT_0323 = IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = IT_0320*IT_0323;
    const ccomplex_t IT_0325 = m_N_4*IT_0229;
    const ccomplex_t IT_0326 = m_N_3*IT_0251;
    const ccomplex_t IT_0327 = IT_0325 + IT_0326;
    const ccomplex_t IT_0328 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0329 = IT_0001*IT_0328;
    const ccomplex_t IT_0330 = 1.4142135623731*IT_0329;
    const ccomplex_t IT_0331 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0332 = IT_0008*IT_0331;
    const ccomplex_t IT_0333 = 1.4142135623731*IT_0332;
    const ccomplex_t IT_0334 = m_t*conj(N_u3)*e_em*IT_0010*U_st_11;
    const ccomplex_t IT_0335 = IT_0009*IT_0334;
    const ccomplex_t IT_0336 = 1.4142135623731*IT_0335;
    const ccomplex_t IT_0337 = (0 + _Complex_I*1)*(IT_0330 + 3*IT_0333 + 3
      *IT_0336);
    const ccomplex_t IT_0338 = 0.166666666666667*IT_0337;
    const ccomplex_t IT_0339 = IT_0206*IT_0338;
    const ccomplex_t IT_0340 = (0 + _Complex_I*1)*IT_0339;
    const ccomplex_t IT_0341 = IT_0327*IT_0340;
    const ccomplex_t IT_0342 = -IT_0269 + -IT_0283 + 2*IT_0287 + -IT_0290 + 
      -IT_0307 + -IT_0324 + -IT_0341;
    const ccomplex_t IT_0343 = 6*s_13;
    const ccomplex_t IT_0344 = m_t*m_N_4;
    const ccomplex_t IT_0345 = 6*IT_0255;
    const ccomplex_t IT_0346 = IT_0067*IT_0289;
    const ccomplex_t IT_0347 = IT_0229*IT_0340;
    const ccomplex_t IT_0348 = IT_0115*IT_0306;
    const ccomplex_t IT_0349 = IT_0194*IT_0320;
    const ccomplex_t IT_0350 = -IT_0053 + IT_0346 + IT_0347 + IT_0348 + IT_0349;
    const ccomplex_t IT_0351 = m_t*s_14;
    const ccomplex_t IT_0352 = IT_0208*IT_0251;
    const ccomplex_t IT_0353 = IT_0077*IT_0082;
    const ccomplex_t IT_0354 = IT_0151*IT_0172;
    const ccomplex_t IT_0355 = IT_0094*IT_0137;
    const ccomplex_t IT_0356 = IT_0286 + -IT_0352 + -IT_0353 + -IT_0354 + 
      -IT_0355;
    const ccomplex_t IT_0357 = s_34*m_N_4;
    const ccomplex_t IT_0358 = 6*conj(IT_0255);
    const ccomplex_t IT_0359 = 6*conj(IT_0342);
    const ccomplex_t IT_0360 = pow(m_H0, 2);
    const ccomplex_t IT_0361 = m_t*m_N_4*IT_0360;
    const ccomplex_t IT_0362 = s_13*IT_0360;
    const ccomplex_t IT_0363 = s_14*s_34;
    const ccomplex_t IT_0364 = (-2)*IT_0363;
    const ccomplex_t IT_0365 = IT_0362 + IT_0364;
    const ccomplex_t IT_0366 = 6*IT_0342;
    const ccomplex_t IT_0367 = 0.166666666666667*IT_0357;
    const ccomplex_t IT_0368 = (IT_0255*conj(IT_0255) + IT_0342*conj(IT_0342))
      *IT_0343 + conj(IT_0342)*IT_0344*IT_0345 + (IT_0342*IT_0344 + IT_0350
      *IT_0351 + IT_0356*IT_0357)*IT_0358 + (IT_0351*IT_0356 + IT_0350*IT_0357)
      *IT_0359 + 6*conj(IT_0356)*(IT_0350*IT_0361 + -IT_0356*IT_0365 +
       0.166666666666667*IT_0351*IT_0366 + IT_0345*IT_0367) + 6*conj(IT_0350)*
      (0.166666666666667*IT_0345*IT_0351 + IT_0356*IT_0361 + -IT_0350*IT_0365 +
       IT_0366*IT_0367);
    return create_ccomplex_return(IT_0368);
}

