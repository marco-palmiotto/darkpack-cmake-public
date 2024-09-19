#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_sb_1_to_s_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_sb_1_to_s_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_ss_R = param->m_ss_R;
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
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_ss_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0004;
    const ccomplex_t IT_0006 = (-0.333333333333333)*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = m_b*conj(N_d2)*e_em*IT_0012*U_sb_10;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0020 = IT_0004*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0015 + -IT_0018 +
       0.333333333333333*IT_0021);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = 3*IT_0024;
    const ccomplex_t IT_0026 = 0.333333333333333*IT_0025;
    const ccomplex_t IT_0027 = 3*IT_0026;
    const ccomplex_t IT_0028 = 0.333333333333333*IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = 3*IT_0032;
    const ccomplex_t IT_0034 = 0.333333333333333*IT_0033;
    const ccomplex_t IT_0035 = 3*IT_0034;
    const ccomplex_t IT_0036 = (-0.333333333333333)*IT_0035;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = 3*IT_0037;
    const ccomplex_t IT_0039 = 0.166666666666667*IT_0038;
    const ccomplex_t IT_0040 = IT_0006*IT_0039;
    const ccomplex_t IT_0041 = IT_0002*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0004;
    const ccomplex_t IT_0045 = (-0.333333333333333)*IT_0044;
    const ccomplex_t IT_0046 = m_b*conj(N_d1)*e_em*IT_0012*U_sb_10;
    const ccomplex_t IT_0047 = IT_0011*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0050 = IT_0004*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0053 = IT_0010*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0048 + 0.333333333333333
      *IT_0051 + -IT_0054);
    const ccomplex_t IT_0056 = 3*IT_0055;
    const ccomplex_t IT_0057 = (-0.333333333333333)*IT_0056;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = 3*IT_0058;
    const ccomplex_t IT_0060 = (-0.333333333333333)*IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = (-3)*IT_0062;
    const ccomplex_t IT_0064 = (-0.333333333333333)*IT_0063;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = 3*IT_0065;
    const ccomplex_t IT_0067 = 0.166666666666667*IT_0066;
    const ccomplex_t IT_0068 = IT_0045*IT_0067;
    const ccomplex_t IT_0069 = IT_0043*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0004;
    const ccomplex_t IT_0073 = (-0.333333333333333)*IT_0072;
    const ccomplex_t IT_0074 = m_b*conj(N_d3)*e_em*IT_0012*U_sb_10;
    const ccomplex_t IT_0075 = IT_0011*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0078 = IT_0010*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0081 = IT_0004*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0076 + -IT_0079 +
       0.333333333333333*IT_0082);
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = 3*IT_0087;
    const ccomplex_t IT_0089 = (-0.333333333333333)*IT_0088;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = 3*IT_0090;
    const ccomplex_t IT_0092 = (-0.333333333333333)*IT_0091;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = 3*IT_0097;
    const ccomplex_t IT_0099 = (-0.333333333333333)*IT_0098;
    const ccomplex_t IT_0100 = (-3)*IT_0099;
    const ccomplex_t IT_0101 = 0.333333333333333*IT_0100;
    const ccomplex_t IT_0102 = 3*IT_0101;
    const ccomplex_t IT_0103 = 0.166666666666667*IT_0102;
    const ccomplex_t IT_0104 = IT_0073*IT_0103;
    const ccomplex_t IT_0105 = IT_0071*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0004;
    const ccomplex_t IT_0109 = (-0.333333333333333)*IT_0108;
    const ccomplex_t IT_0110 = m_b*conj(N_d4)*e_em*IT_0012*U_sb_10;
    const ccomplex_t IT_0111 = IT_0011*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0114 = IT_0010*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0117 = IT_0004*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0112 + -IT_0115 +
       0.333333333333333*IT_0118);
    const ccomplex_t IT_0120 = 3*IT_0119;
    const ccomplex_t IT_0121 = 0.333333333333333*IT_0120;
    const ccomplex_t IT_0122 = 3*IT_0121;
    const ccomplex_t IT_0123 = 0.166666666666667*IT_0122;
    const ccomplex_t IT_0124 = IT_0109*IT_0123;
    const ccomplex_t IT_0125 = IT_0107*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = pow(m_sG, 2);
    const ccomplex_t IT_0128 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0127 +
       reg_prop, -1);
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0131 = IT_0129*IT_0130;
    const ccomplex_t IT_0132 = IT_0128*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = -IT_0042 + -IT_0070 + -IT_0106 + -IT_0126 +
       0.166666666666667*IT_0133;
    const ccomplex_t IT_0135 = (-0.5)*IT_0133;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0139 = IT_0004*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = m_b*N_d1*e_em*IT_0012*U_sb_00;
    const ccomplex_t IT_0142 = IT_0011*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0140 + 1.5*IT_0143);
    const ccomplex_t IT_0145 = (-0.333333333333333)*IT_0144;
    const ccomplex_t IT_0146 = IT_0137*IT_0145;
    const ccomplex_t IT_0147 = IT_0043*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0152 = IT_0004*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = m_b*N_d2*e_em*IT_0012*U_sb_00;
    const ccomplex_t IT_0155 = IT_0011*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0153 + 1.5*IT_0156);
    const ccomplex_t IT_0158 = 0.666666666666667*IT_0157;
    const ccomplex_t IT_0159 = 1.5*IT_0158;
    const ccomplex_t IT_0160 = 0.666666666666667*IT_0159;
    const ccomplex_t IT_0161 = 1.5*IT_0160;
    const ccomplex_t IT_0162 = 0.666666666666667*IT_0161;
    const ccomplex_t IT_0163 = 1.5*IT_0162;
    const ccomplex_t IT_0164 = 0.666666666666667*IT_0163;
    const ccomplex_t IT_0165 = 1.5*IT_0164;
    const ccomplex_t IT_0166 = 0.666666666666667*IT_0165;
    const ccomplex_t IT_0167 = 1.5*IT_0166;
    const ccomplex_t IT_0168 = 0.666666666666667*IT_0167;
    const ccomplex_t IT_0169 = 1.5*IT_0168;
    const ccomplex_t IT_0170 = 0.666666666666667*IT_0169;
    const ccomplex_t IT_0171 = 1.5*IT_0170;
    const ccomplex_t IT_0172 = (-0.333333333333333)*IT_0171;
    const ccomplex_t IT_0173 = IT_0150*IT_0172;
    const ccomplex_t IT_0174 = IT_0002*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0177 = 0.5*IT_0176;
    const ccomplex_t IT_0178 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0179 = IT_0004*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = m_b*N_d3*e_em*IT_0012*U_sb_00;
    const ccomplex_t IT_0182 = IT_0011*IT_0181;
    const ccomplex_t IT_0183 = 1.4142135623731*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*(IT_0180 + 1.5*IT_0183);
    const ccomplex_t IT_0185 = 0.666666666666667*IT_0184;
    const ccomplex_t IT_0186 = 1.5*IT_0185;
    const ccomplex_t IT_0187 = (-0.333333333333333)*IT_0186;
    const ccomplex_t IT_0188 = IT_0177*IT_0187;
    const ccomplex_t IT_0189 = IT_0071*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0194 = IT_0004*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = m_b*N_d4*e_em*IT_0012*U_sb_00;
    const ccomplex_t IT_0197 = IT_0011*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0195 + 1.5*IT_0198);
    const ccomplex_t IT_0200 = 0.666666666666667*IT_0199;
    const ccomplex_t IT_0201 = 1.5*IT_0200;
    const ccomplex_t IT_0202 = (-0.333333333333333)*IT_0201;
    const ccomplex_t IT_0203 = IT_0192*IT_0202;
    const ccomplex_t IT_0204 = IT_0107*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = -IT_0148 + -IT_0175 + -IT_0190 + -IT_0205;
    const ccomplex_t IT_0207 = m_b*m_s*IT_0001;
    const ccomplex_t IT_0208 = s_34*IT_0001;
    const ccomplex_t IT_0209 = (-18)*IT_0208;
    const ccomplex_t IT_0210 = s_13*s_14;
    const ccomplex_t IT_0211 = 36*IT_0210;
    const ccomplex_t IT_0212 = IT_0209 + IT_0211;
    const ccomplex_t IT_0213 = m_s*s_14;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0215 = IT_0067*IT_0137;
    const ccomplex_t IT_0216 = IT_0043*IT_0214*IT_0215;
    const ccomplex_t IT_0217 = m_s*IT_0070;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0219 = IT_0123*IT_0192;
    const ccomplex_t IT_0220 = IT_0107*IT_0218*IT_0219;
    const ccomplex_t IT_0221 = m_s*IT_0133;
    const ccomplex_t IT_0222 = m_s*IT_0042;
    const ccomplex_t IT_0223 = m_s*IT_0106;
    const ccomplex_t IT_0224 = m_s*IT_0126;
    const ccomplex_t IT_0225 = IT_0216 + -IT_0217 + IT_0220 +
       0.166666666666667*IT_0221 + -IT_0222 + -IT_0223 + -IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0227 = IT_0039*IT_0150;
    const ccomplex_t IT_0228 = IT_0002*IT_0226*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0230 = IT_0103*IT_0177;
    const ccomplex_t IT_0231 = IT_0071*IT_0229*IT_0230;
    const ccomplex_t IT_0232 = -IT_0228 + -IT_0231;
    const ccomplex_t IT_0233 = (-0.5)*IT_0221;
    const ccomplex_t IT_0234 = -IT_0232;
    const ccomplex_t IT_0235 = IT_0225 + IT_0234;
    const ccomplex_t IT_0236 = -conj(IT_0232);
    const ccomplex_t IT_0237 = conj(IT_0225) + IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0239 = IT_0129*IT_0238;
    const ccomplex_t IT_0240 = IT_0128*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = m_sG*IT_0241;
    const ccomplex_t IT_0243 = (-0.5)*IT_0242;
    const ccomplex_t IT_0244 = 6*IT_0243;
    const ccomplex_t IT_0245 = IT_0045*IT_0145;
    const ccomplex_t IT_0246 = IT_0043*IT_0214*IT_0245;
    const ccomplex_t IT_0247 = 0.166666666666667*IT_0242 + -IT_0246;
    const ccomplex_t IT_0248 = 18*IT_0247;
    const ccomplex_t IT_0249 = IT_0073*IT_0187;
    const ccomplex_t IT_0250 = IT_0071*IT_0229*IT_0249;
    const ccomplex_t IT_0251 = IT_0109*IT_0202;
    const ccomplex_t IT_0252 = IT_0107*IT_0218*IT_0251;
    const ccomplex_t IT_0253 = -IT_0250 + -IT_0252;
    const ccomplex_t IT_0254 = 18*IT_0253;
    const ccomplex_t IT_0255 = IT_0006*IT_0172;
    const ccomplex_t IT_0256 = IT_0002*IT_0226*IT_0255;
    const ccomplex_t IT_0257 = m_s*IT_0148;
    const ccomplex_t IT_0258 = m_s*IT_0175;
    const ccomplex_t IT_0259 = m_s*IT_0190;
    const ccomplex_t IT_0260 = m_s*IT_0205;
    const ccomplex_t IT_0261 = IT_0256 + -IT_0257 + -IT_0258 + -IT_0259 + 
      -IT_0260;
    const ccomplex_t IT_0262 = (-18)*IT_0261;
    const ccomplex_t IT_0263 = IT_0244 + IT_0248 + IT_0254 + IT_0262;
    const ccomplex_t IT_0264 = 6*conj(IT_0243);
    const ccomplex_t IT_0265 = 18*conj(IT_0247);
    const ccomplex_t IT_0266 = 18*conj(IT_0253);
    const ccomplex_t IT_0267 = (-18)*conj(IT_0261);
    const ccomplex_t IT_0268 = IT_0264 + IT_0265 + IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = conj(IT_0134)*((-18)*IT_0225 + 18*IT_0232 + (-6
      )*IT_0233) + IT_0134*((-18)*conj(IT_0225) + 18*conj(IT_0232) + (-6)*conj
      (IT_0233)) + (-18)*conj(IT_0135)*(IT_0233 + 0.333333333333333*IT_0235) + (
      -18)*IT_0135*(conj(IT_0233) + 0.333333333333333*IT_0237) + conj(IT_0206)
      *IT_0263 + IT_0206*IT_0268;
    const ccomplex_t IT_0270 = conj(IT_0247) + conj(IT_0253);
    const ccomplex_t IT_0271 = 6*s_34;
    const ccomplex_t IT_0272 = 18*s_34;
    const ccomplex_t IT_0273 = (-18)*s_34;
    const ccomplex_t IT_0274 = (-6)*s_34;
    const ccomplex_t IT_0275 = 12*IT_0210;
    const ccomplex_t IT_0276 = (-6)*IT_0208;
    const ccomplex_t IT_0277 = m_b*m_s;
    const ccomplex_t IT_0278 = 18*IT_0261;
    const ccomplex_t IT_0279 = (-18)*IT_0247;
    const ccomplex_t IT_0280 = (-18)*IT_0253;
    const ccomplex_t IT_0281 = (-6)*IT_0243;
    const ccomplex_t IT_0282 = IT_0278 + IT_0279 + IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = 18*conj(IT_0261);
    const ccomplex_t IT_0284 = (-18)*conj(IT_0247);
    const ccomplex_t IT_0285 = (-18)*conj(IT_0253);
    const ccomplex_t IT_0286 = (-6)*conj(IT_0243);
    const ccomplex_t IT_0287 = IT_0283 + IT_0284 + IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = 3*conj(IT_0243);
    const ccomplex_t IT_0289 = -conj(IT_0261);
    const ccomplex_t IT_0290 = 3*conj(IT_0233)*(IT_0243 + 1./3*IT_0247 + 1./3
      *IT_0253 + -1./3*IT_0261) + 0.166666666666667*conj(IT_0225)*IT_0263 +
       0.166666666666667*IT_0225*IT_0268 + 0.166666666666667*conj(IT_0232)
      *IT_0282 + 0.166666666666667*IT_0232*IT_0287 + IT_0233*(conj(IT_0247) +
       conj(IT_0253) + IT_0288 + IT_0289);
    const ccomplex_t IT_0291 = m_b*s_13;
    const ccomplex_t IT_0292 = -conj(IT_0206)*(IT_0233 + 3*IT_0235) + -IT_0206
      *(conj(IT_0233) + 3*IT_0237) + 3*conj(IT_0135)*(IT_0243 + 1./3*IT_0247 +
       1./3*IT_0253 + -1./3*IT_0261) + (-0.166666666666667)*conj(IT_0134)
      *IT_0282 + (-0.166666666666667)*IT_0134*IT_0287 + IT_0135*(conj(IT_0247) +
       conj(IT_0253) + IT_0288 + IT_0289);
    const ccomplex_t IT_0293 = (((-18)*conj(IT_0134) + (-6)*conj(IT_0135))
      *IT_0206 + ((-18)*IT_0134 + (-6)*IT_0135)*conj(IT_0206))*IT_0207 + 
      (IT_0134*conj(IT_0134) + IT_0135*conj(IT_0135) + IT_0206*conj(IT_0206))
      *IT_0212 + IT_0213*IT_0269 + (conj(IT_0225)*IT_0233 + IT_0225*conj(IT_0233
      ) + conj(IT_0243)*(IT_0247 + IT_0253) + IT_0243*IT_0270)*IT_0271 + 
      (IT_0225*conj(IT_0225) + IT_0232*conj(IT_0232) + IT_0233*conj(IT_0233) +
       IT_0243*conj(IT_0243) + IT_0261*conj(IT_0261) + (IT_0247 + IT_0253)
      *IT_0270)*IT_0272 + (conj(IT_0225)*IT_0232 + IT_0225*conj(IT_0232) + 
      (IT_0247 + IT_0253)*conj(IT_0261) + IT_0261*IT_0270)*IT_0273 + (conj
      (IT_0232)*IT_0233 + IT_0232*conj(IT_0233) + conj(IT_0243)*IT_0261 +
       IT_0243*conj(IT_0261))*IT_0274 + (conj(IT_0134)*IT_0135 + IT_0134*conj
      (IT_0135))*(IT_0275 + IT_0276) + 6*IT_0277*IT_0290 + (-6)*IT_0291*IT_0292;
    return create_ccomplex_return(IT_0293);
}

