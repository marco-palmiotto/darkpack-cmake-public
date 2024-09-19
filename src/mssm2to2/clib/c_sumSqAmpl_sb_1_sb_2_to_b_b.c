#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_sb_2_to_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_sb_2_to_b_b(
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
    const creal_t m_sb_2 = param->m_sb_2;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_sb_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
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
    const ccomplex_t IT_0025 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0026 = IT_0010*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0029 = IT_0017*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = m_b*conj(N_d1)*e_em*IT_0019*U_sb_11;
    const ccomplex_t IT_0032 = IT_0018*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0027 + (-3)*IT_0030 + 3
      *IT_0033);
    const ccomplex_t IT_0035 = 0.166666666666667*IT_0034;
    const ccomplex_t IT_0036 = IT_0024*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0041 = IT_0010*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = m_b*N_d2*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0044 = IT_0018*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 1.5*IT_0045);
    const ccomplex_t IT_0047 = (-0.333333333333333)*IT_0046;
    const ccomplex_t IT_0048 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0049 = IT_0010*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0052 = IT_0017*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = m_b*conj(N_d2)*e_em*IT_0019*U_sb_11;
    const ccomplex_t IT_0055 = IT_0018*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0050 + (-3)*IT_0053 + 3
      *IT_0056);
    const ccomplex_t IT_0058 = 0.166666666666667*IT_0057;
    const ccomplex_t IT_0059 = IT_0047*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0064 = IT_0010*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = m_b*N_d3*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0067 = IT_0018*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0065 + 1.5*IT_0068);
    const ccomplex_t IT_0070 = (-0.333333333333333)*IT_0069;
    const ccomplex_t IT_0071 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0072 = IT_0010*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0075 = IT_0017*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = m_b*conj(N_d3)*e_em*IT_0019*U_sb_11;
    const ccomplex_t IT_0078 = IT_0018*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0073 + (-3)*IT_0076 + 3
      *IT_0079);
    const ccomplex_t IT_0081 = 0.166666666666667*IT_0080;
    const ccomplex_t IT_0082 = IT_0070*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0087 = IT_0010*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = m_b*N_d4*e_em*IT_0019*U_sb_00;
    const ccomplex_t IT_0090 = IT_0018*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + 1.5*IT_0091);
    const ccomplex_t IT_0093 = (-0.333333333333333)*IT_0092;
    const ccomplex_t IT_0094 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0095 = IT_0010*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0098 = IT_0017*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = m_b*conj(N_d4)*e_em*IT_0019*U_sb_11;
    const ccomplex_t IT_0101 = IT_0018*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0096 + (-3)*IT_0099 + 3
      *IT_0102);
    const ccomplex_t IT_0104 = 0.166666666666667*IT_0103;
    const ccomplex_t IT_0105 = IT_0093*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0107 = IT_0105*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = pow(m_sb_1, 2);
    const ccomplex_t IT_0110 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0109 +
       reg_prop, -1);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0113 = IT_0111*IT_0112;
    const ccomplex_t IT_0114 = IT_0110*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = (-0.166666666666667)*IT_0008 + IT_0039 +
       IT_0062 + IT_0085 + IT_0108 + (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0118 = IT_0010*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0121 = IT_0017*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = m_b*conj(N_d1)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0124 = IT_0018*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0119 + (-3)*IT_0122 + 3
      *IT_0125);
    const ccomplex_t IT_0127 = 0.166666666666667*IT_0126;
    const ccomplex_t IT_0128 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0129 = IT_0010*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = m_b*N_d1*e_em*IT_0019*U_sb_01;
    const ccomplex_t IT_0132 = IT_0018*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0130 + 1.5*IT_0133);
    const ccomplex_t IT_0135 = (-0.333333333333333)*IT_0134;
    const ccomplex_t IT_0136 = IT_0127*IT_0135;
    const ccomplex_t IT_0137 = cpow((-2)*s_13 + IT_0000 + IT_0109 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0138 = IT_0136*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0141 = IT_0010*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0144 = IT_0017*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = m_b*conj(N_d2)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0147 = IT_0018*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0142 + (-3)*IT_0145 + 3
      *IT_0148);
    const ccomplex_t IT_0150 = 0.166666666666667*IT_0149;
    const ccomplex_t IT_0151 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0152 = IT_0010*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = m_b*N_d2*e_em*IT_0019*U_sb_01;
    const ccomplex_t IT_0155 = IT_0018*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0153 + 1.5*IT_0156);
    const ccomplex_t IT_0158 = (-0.333333333333333)*IT_0157;
    const ccomplex_t IT_0159 = IT_0150*IT_0158;
    const ccomplex_t IT_0160 = cpow((-2)*s_13 + IT_0000 + IT_0109 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0161 = IT_0159*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0164 = IT_0010*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0167 = IT_0017*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = m_b*conj(N_d3)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0170 = IT_0018*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*(IT_0165 + (-3)*IT_0168 + 3
      *IT_0171);
    const ccomplex_t IT_0173 = 0.166666666666667*IT_0172;
    const ccomplex_t IT_0174 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0175 = IT_0010*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = m_b*N_d3*e_em*IT_0019*U_sb_01;
    const ccomplex_t IT_0178 = IT_0018*IT_0177;
    const ccomplex_t IT_0179 = 1.4142135623731*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*(IT_0176 + 1.5*IT_0179);
    const ccomplex_t IT_0181 = (-0.333333333333333)*IT_0180;
    const ccomplex_t IT_0182 = IT_0173*IT_0181;
    const ccomplex_t IT_0183 = cpow((-2)*s_13 + IT_0000 + IT_0109 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0184 = IT_0182*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0187 = IT_0010*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0190 = IT_0017*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = m_b*conj(N_d4)*e_em*IT_0019*U_sb_10;
    const ccomplex_t IT_0193 = IT_0018*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*(IT_0188 + (-3)*IT_0191 + 3
      *IT_0194);
    const ccomplex_t IT_0196 = 0.166666666666667*IT_0195;
    const ccomplex_t IT_0197 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0198 = IT_0010*IT_0197;
    const ccomplex_t IT_0199 = 1.4142135623731*IT_0198;
    const ccomplex_t IT_0200 = m_b*N_d4*e_em*IT_0019*U_sb_01;
    const ccomplex_t IT_0201 = IT_0018*IT_0200;
    const ccomplex_t IT_0202 = 1.4142135623731*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*(IT_0199 + 1.5*IT_0202);
    const ccomplex_t IT_0204 = (-0.333333333333333)*IT_0203;
    const ccomplex_t IT_0205 = IT_0196*IT_0204;
    const ccomplex_t IT_0206 = cpow((-2)*s_13 + IT_0000 + IT_0109 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0207 = IT_0205*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = 0.5*IT_0008 + 0.166666666666667*IT_0115 + 
      -IT_0139 + -IT_0162 + -IT_0185 + -IT_0208;
    const ccomplex_t IT_0210 = IT_0006*IT_0110;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = IT_0003*IT_0113;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = IT_0037*IT_0136;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = IT_0060*IT_0159;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = IT_0083*IT_0182;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = IT_0106*IT_0205;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = (-0.5)*IT_0211 + (-0.166666666666667)*IT_0213 +
       IT_0215 + IT_0217 + IT_0219 + IT_0221;
    const ccomplex_t IT_0223 = IT_0036*IT_0137;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = IT_0059*IT_0160;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = IT_0082*IT_0183;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = IT_0105*IT_0206;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = 0.166666666666667*IT_0211 + 0.5*IT_0213 + 
      -IT_0224 + -IT_0226 + -IT_0228 + -IT_0230;
    const ccomplex_t IT_0232 = s_34*IT_0002;
    const ccomplex_t IT_0233 = (-6)*IT_0232;
    const ccomplex_t IT_0234 = s_23*s_24;
    const ccomplex_t IT_0235 = 12*IT_0234;
    const ccomplex_t IT_0236 = IT_0233 + IT_0235;
    const ccomplex_t IT_0237 = (-18)*IT_0232;
    const ccomplex_t IT_0238 = 36*IT_0234;
    const ccomplex_t IT_0239 = IT_0237 + IT_0238;
    const ccomplex_t IT_0240 = m_b*s_24;
    const ccomplex_t IT_0241 = IT_0004*IT_0112;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0243 = IT_0110*IT_0241*IT_0242;
    const ccomplex_t IT_0244 = IT_0003*IT_0241*IT_0242;
    const ccomplex_t IT_0245 = IT_0024*IT_0135;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0247 = IT_0037*IT_0245*IT_0246;
    const ccomplex_t IT_0248 = IT_0047*IT_0158;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0250 = IT_0060*IT_0248*IT_0249;
    const ccomplex_t IT_0251 = IT_0070*IT_0181;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0253 = IT_0083*IT_0251*IT_0252;
    const ccomplex_t IT_0254 = IT_0093*IT_0204;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0256 = IT_0106*IT_0254*IT_0255;
    const ccomplex_t IT_0257 = m_b*IT_0211;
    const ccomplex_t IT_0258 = m_b*IT_0008;
    const ccomplex_t IT_0259 = m_b*IT_0039;
    const ccomplex_t IT_0260 = m_b*IT_0062;
    const ccomplex_t IT_0261 = m_b*IT_0085;
    const ccomplex_t IT_0262 = m_b*IT_0108;
    const ccomplex_t IT_0263 = 0.5*IT_0243 + (-0.166666666666667)*IT_0244 +
       IT_0247 + IT_0250 + IT_0253 + IT_0256 + (-0.5)*IT_0257 +
       0.166666666666667*IT_0258 + -IT_0259 + -IT_0260 + -IT_0261 + -IT_0262;
    const ccomplex_t IT_0264 = IT_0005*IT_0111;
    const ccomplex_t IT_0265 = IT_0110*IT_0242*IT_0264;
    const ccomplex_t IT_0266 = IT_0003*IT_0242*IT_0264;
    const ccomplex_t IT_0267 = IT_0035*IT_0127;
    const ccomplex_t IT_0268 = IT_0037*IT_0246*IT_0267;
    const ccomplex_t IT_0269 = IT_0058*IT_0150;
    const ccomplex_t IT_0270 = IT_0060*IT_0249*IT_0269;
    const ccomplex_t IT_0271 = IT_0081*IT_0173;
    const ccomplex_t IT_0272 = IT_0083*IT_0252*IT_0271;
    const ccomplex_t IT_0273 = IT_0104*IT_0196;
    const ccomplex_t IT_0274 = IT_0106*IT_0255*IT_0273;
    const ccomplex_t IT_0275 = m_b*IT_0213;
    const ccomplex_t IT_0276 = m_b*IT_0115;
    const ccomplex_t IT_0277 = m_b*IT_0215;
    const ccomplex_t IT_0278 = m_b*IT_0217;
    const ccomplex_t IT_0279 = m_b*IT_0219;
    const ccomplex_t IT_0280 = m_b*IT_0221;
    const ccomplex_t IT_0281 = 0.5*IT_0265 + (-0.166666666666667)*IT_0266 +
       IT_0268 + IT_0270 + IT_0272 + IT_0274 + 0.166666666666667*IT_0275 + (-0.5
      )*IT_0276 + -IT_0277 + -IT_0278 + -IT_0279 + -IT_0280;
    const ccomplex_t IT_0282 = IT_0137*IT_0245*IT_0246;
    const ccomplex_t IT_0283 = IT_0160*IT_0248*IT_0249;
    const ccomplex_t IT_0284 = IT_0183*IT_0251*IT_0252;
    const ccomplex_t IT_0285 = IT_0206*IT_0254*IT_0255;
    const ccomplex_t IT_0286 = m_b*IT_0224;
    const ccomplex_t IT_0287 = m_b*IT_0226;
    const ccomplex_t IT_0288 = m_b*IT_0228;
    const ccomplex_t IT_0289 = m_b*IT_0230;
    const ccomplex_t IT_0290 = (-0.166666666666667)*IT_0243 + 0.5*IT_0244 +
       0.166666666666667*IT_0257 + (-0.5)*IT_0258 + IT_0282 + IT_0283 + IT_0284 
      + IT_0285 + -IT_0286 + -IT_0287 + -IT_0288 + -IT_0289;
    const ccomplex_t IT_0291 = IT_0116 + 3*IT_0209;
    const ccomplex_t IT_0292 = 3*conj(IT_0209);
    const ccomplex_t IT_0293 = conj(IT_0116) + IT_0292;
    const ccomplex_t IT_0294 = IT_0137*IT_0246*IT_0267;
    const ccomplex_t IT_0295 = IT_0160*IT_0249*IT_0269;
    const ccomplex_t IT_0296 = IT_0183*IT_0252*IT_0271;
    const ccomplex_t IT_0297 = IT_0206*IT_0255*IT_0273;
    const ccomplex_t IT_0298 = m_b*IT_0139;
    const ccomplex_t IT_0299 = m_b*IT_0162;
    const ccomplex_t IT_0300 = m_b*IT_0185;
    const ccomplex_t IT_0301 = m_b*IT_0208;
    const ccomplex_t IT_0302 = (-0.166666666666667)*IT_0265 + 0.5*IT_0266 + (
      -0.5)*IT_0275 + 0.166666666666667*IT_0276 + IT_0294 + IT_0295 + IT_0296 +
       IT_0297 + -IT_0298 + -IT_0299 + -IT_0300 + -IT_0301;
    const ccomplex_t IT_0303 = 3*IT_0231;
    const ccomplex_t IT_0304 = IT_0222 + IT_0303;
    const ccomplex_t IT_0305 = 3*conj(IT_0231);
    const ccomplex_t IT_0306 = conj(IT_0222) + IT_0305;
    const ccomplex_t IT_0307 = (18*conj(IT_0116) + 6*conj(IT_0209))*IT_0263 + 
      (18*IT_0116 + 6*IT_0209)*conj(IT_0263) + (18*conj(IT_0222) + 6*conj
      (IT_0231))*IT_0281 + (18*IT_0222 + 6*IT_0231)*conj(IT_0281) + 6*conj
      (IT_0290)*IT_0291 + 6*IT_0290*IT_0293 + 6*conj(IT_0302)*IT_0304 + 6
      *IT_0302*IT_0306;
    const ccomplex_t IT_0308 = 6*s_34;
    const ccomplex_t IT_0309 = 18*s_34;
    const ccomplex_t IT_0310 = (-18)*conj(IT_0263);
    const ccomplex_t IT_0311 = (-6)*conj(IT_0290);
    const ccomplex_t IT_0312 = IT_0310 + IT_0311;
    const ccomplex_t IT_0313 = (-18)*conj(IT_0290);
    const ccomplex_t IT_0314 = (-6)*conj(IT_0263);
    const ccomplex_t IT_0315 = IT_0313 + IT_0314;
    const ccomplex_t IT_0316 = (-18)*IT_0290;
    const ccomplex_t IT_0317 = (-6)*IT_0290;
    const ccomplex_t IT_0318 = IT_0000*IT_0002;
    const ccomplex_t IT_0319 = (-6)*conj(IT_0231);
    const ccomplex_t IT_0320 = (-18)*conj(IT_0222);
    const ccomplex_t IT_0321 = IT_0319 + IT_0320;
    const ccomplex_t IT_0322 = (-6)*IT_0209;
    const ccomplex_t IT_0323 = (-6)*conj(IT_0209);
    const ccomplex_t IT_0324 = m_b*s_23;
    const ccomplex_t IT_0325 = IT_0293*IT_0302 + IT_0291*conj(IT_0302) + (
      -0.166666666666667)*IT_0222*IT_0312 + (-0.166666666666667)*IT_0231*IT_0315
       + (-0.166666666666667)*conj(IT_0231)*IT_0316 + (-0.166666666666667)*conj
      (IT_0222)*IT_0317 + (-0.166666666666667)*IT_0263*IT_0321 + (
      -0.166666666666667)*conj(IT_0281)*((-18)*IT_0116 + IT_0322) + (
      -0.166666666666667)*IT_0281*((-18)*conj(IT_0116) + IT_0323);
    const ccomplex_t IT_0326 = (conj(IT_0116)*IT_0209 + IT_0116*conj(IT_0209) 
      + conj(IT_0222)*IT_0231 + IT_0222*conj(IT_0231))*IT_0236 + (IT_0116*conj
      (IT_0116) + IT_0209*conj(IT_0209) + IT_0222*conj(IT_0222) + IT_0231*conj
      (IT_0231))*IT_0239 + IT_0240*IT_0307 + (conj(IT_0263)*IT_0290 + IT_0263
      *conj(IT_0290) + conj(IT_0281)*IT_0302 + IT_0281*conj(IT_0302))*IT_0308 + 
      (IT_0263*conj(IT_0263) + IT_0281*conj(IT_0281) + IT_0290*conj(IT_0290) +
       IT_0302*conj(IT_0302))*IT_0309 + IT_0000*(IT_0263*((-18)*conj(IT_0281) + 
      (-6)*conj(IT_0302)) + IT_0281*IT_0312 + IT_0302*IT_0315 + conj(IT_0302)
      *IT_0316 + conj(IT_0281)*IT_0317) + IT_0318*(conj(IT_0116)*((-18)*IT_0222 
      + (-6)*IT_0231) + IT_0116*IT_0321 + IT_0306*IT_0322 + IT_0304*IT_0323) + (
      -6)*IT_0324*IT_0325;
    return create_ccomplex_return(IT_0326);
}

