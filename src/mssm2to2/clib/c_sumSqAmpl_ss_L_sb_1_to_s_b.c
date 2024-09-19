#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_sb_1_to_s_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_sb_1_to_s_b(
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
    const creal_t m_ss_L = param->m_ss_L;
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
    const ccomplex_t IT_0001 = pow(m_ss_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = IT_0005*IT_0007;
    const ccomplex_t IT_0011 = m_b*conj(N_d1)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = cos(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0013 + -IT_0016 +
       0.333333333333333*IT_0021);
    const ccomplex_t IT_0023 = 3*IT_0022;
    const ccomplex_t IT_0024 = (-0.333333333333333)*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = 3*IT_0025;
    const ccomplex_t IT_0027 = (-0.333333333333333)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = (-3)*IT_0029;
    const ccomplex_t IT_0031 = (-0.333333333333333)*IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = 3*IT_0032;
    const ccomplex_t IT_0034 = 0.166666666666667*IT_0033;
    const ccomplex_t IT_0035 = IT_0009*IT_0034;
    const ccomplex_t IT_0036 = IT_0002*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = m_b*conj(N_d2)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0042 = IT_0010*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0045 = IT_0007*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0048 = IT_0018*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0043 + -IT_0046 +
       0.333333333333333*IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = 3*IT_0052;
    const ccomplex_t IT_0054 = 0.333333333333333*IT_0053;
    const ccomplex_t IT_0055 = 3*IT_0054;
    const ccomplex_t IT_0056 = 0.333333333333333*IT_0055;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = 3*IT_0060;
    const ccomplex_t IT_0062 = 0.333333333333333*IT_0061;
    const ccomplex_t IT_0063 = 3*IT_0062;
    const ccomplex_t IT_0064 = (-0.333333333333333)*IT_0063;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = 3*IT_0065;
    const ccomplex_t IT_0067 = 0.166666666666667*IT_0066;
    const ccomplex_t IT_0068 = IT_0040*IT_0067;
    const ccomplex_t IT_0069 = IT_0038*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = m_b*conj(N_d3)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0075 = IT_0010*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0078 = IT_0018*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0081 = IT_0007*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0076 + 0.333333333333333
      *IT_0079 + -IT_0082);
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
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = m_b*conj(N_d4)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0111 = IT_0010*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0114 = IT_0007*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0117 = IT_0018*IT_0116;
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
    const ccomplex_t IT_0127 = -IT_0037 + -IT_0070 + -IT_0106 + -IT_0126;
    const ccomplex_t IT_0128 = conj(N_B1)*e_em;
    const ccomplex_t IT_0129 = IT_0018*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = conj(N_W1)*e_em;
    const ccomplex_t IT_0132 = IT_0007*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0130 + (-3)*IT_0133);
    const ccomplex_t IT_0135 = 0.166666666666667*IT_0134;
    const ccomplex_t IT_0136 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0137 = IT_0018*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = m_b*N_d1*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0140 = IT_0010*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0138 + 1.5*IT_0141);
    const ccomplex_t IT_0143 = (-0.333333333333333)*IT_0142;
    const ccomplex_t IT_0144 = IT_0135*IT_0143;
    const ccomplex_t IT_0145 = IT_0002*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = pow(m_sG, 2);
    const ccomplex_t IT_0148 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0147 +
       reg_prop, -1);
    const ccomplex_t IT_0149 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0151 = IT_0149*IT_0150;
    const ccomplex_t IT_0152 = IT_0148*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = conj(N_B2)*e_em;
    const ccomplex_t IT_0155 = IT_0018*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = conj(N_W2)*e_em;
    const ccomplex_t IT_0158 = IT_0007*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*(IT_0156 + (-3)*IT_0159);
    const ccomplex_t IT_0161 = 0.166666666666667*IT_0160;
    const ccomplex_t IT_0162 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0163 = IT_0018*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = m_b*N_d2*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0166 = IT_0010*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*(IT_0164 + 1.5*IT_0167);
    const ccomplex_t IT_0169 = 0.666666666666667*IT_0168;
    const ccomplex_t IT_0170 = 1.5*IT_0169;
    const ccomplex_t IT_0171 = 0.666666666666667*IT_0170;
    const ccomplex_t IT_0172 = 1.5*IT_0171;
    const ccomplex_t IT_0173 = 0.666666666666667*IT_0172;
    const ccomplex_t IT_0174 = 1.5*IT_0173;
    const ccomplex_t IT_0175 = 0.666666666666667*IT_0174;
    const ccomplex_t IT_0176 = 1.5*IT_0175;
    const ccomplex_t IT_0177 = 0.666666666666667*IT_0176;
    const ccomplex_t IT_0178 = 1.5*IT_0177;
    const ccomplex_t IT_0179 = 0.666666666666667*IT_0178;
    const ccomplex_t IT_0180 = 1.5*IT_0179;
    const ccomplex_t IT_0181 = 0.666666666666667*IT_0180;
    const ccomplex_t IT_0182 = 1.5*IT_0181;
    const ccomplex_t IT_0183 = (-0.333333333333333)*IT_0182;
    const ccomplex_t IT_0184 = IT_0161*IT_0183;
    const ccomplex_t IT_0185 = IT_0038*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = conj(N_B3)*e_em;
    const ccomplex_t IT_0188 = IT_0018*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = conj(N_W3)*e_em;
    const ccomplex_t IT_0191 = IT_0007*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*(IT_0189 + (-3)*IT_0192);
    const ccomplex_t IT_0194 = 0.166666666666667*IT_0193;
    const ccomplex_t IT_0195 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0196 = IT_0018*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = m_b*N_d3*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0199 = IT_0010*IT_0198;
    const ccomplex_t IT_0200 = 1.4142135623731*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*(IT_0197 + 1.5*IT_0200);
    const ccomplex_t IT_0202 = 0.666666666666667*IT_0201;
    const ccomplex_t IT_0203 = 1.5*IT_0202;
    const ccomplex_t IT_0204 = (-0.333333333333333)*IT_0203;
    const ccomplex_t IT_0205 = IT_0194*IT_0204;
    const ccomplex_t IT_0206 = IT_0071*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*IT_0206;
    const ccomplex_t IT_0208 = conj(N_B4)*e_em;
    const ccomplex_t IT_0209 = IT_0018*IT_0208;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = conj(N_W4)*e_em;
    const ccomplex_t IT_0212 = IT_0007*IT_0211;
    const ccomplex_t IT_0213 = 1.4142135623731*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*(IT_0210 + (-3)*IT_0213);
    const ccomplex_t IT_0215 = 0.166666666666667*IT_0214;
    const ccomplex_t IT_0216 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0217 = IT_0018*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = m_b*N_d4*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0220 = IT_0010*IT_0219;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*(IT_0218 + 1.5*IT_0221);
    const ccomplex_t IT_0223 = 0.666666666666667*IT_0222;
    const ccomplex_t IT_0224 = 1.5*IT_0223;
    const ccomplex_t IT_0225 = (-0.333333333333333)*IT_0224;
    const ccomplex_t IT_0226 = IT_0215*IT_0225;
    const ccomplex_t IT_0227 = IT_0107*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = -IT_0146 + 0.166666666666667*IT_0153 + -IT_0186
       + -IT_0207 + -IT_0228;
    const ccomplex_t IT_0230 = (-0.5)*IT_0153;
    const ccomplex_t IT_0231 = m_b*m_s*IT_0001;
    const ccomplex_t IT_0232 = s_34*IT_0001;
    const ccomplex_t IT_0233 = (-18)*IT_0232;
    const ccomplex_t IT_0234 = s_13*s_14;
    const ccomplex_t IT_0235 = 36*IT_0234;
    const ccomplex_t IT_0236 = IT_0233 + IT_0235;
    const ccomplex_t IT_0237 = m_b*s_13;
    const ccomplex_t IT_0238 = m_s*IT_0126;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0240 = IT_0103*IT_0194;
    const ccomplex_t IT_0241 = IT_0071*IT_0239*IT_0240;
    const ccomplex_t IT_0242 = m_s*IT_0037;
    const ccomplex_t IT_0243 = m_s*IT_0070;
    const ccomplex_t IT_0244 = m_s*IT_0106;
    const ccomplex_t IT_0245 = -IT_0238 + IT_0241 + -IT_0242 + -IT_0243 + 
      -IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0247 = IT_0149*IT_0246;
    const ccomplex_t IT_0248 = IT_0148*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = m_sG*IT_0249;
    const ccomplex_t IT_0251 = (-0.5)*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0253 = IT_0034*IT_0135;
    const ccomplex_t IT_0254 = IT_0002*IT_0252*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0256 = IT_0067*IT_0161;
    const ccomplex_t IT_0257 = IT_0038*IT_0255*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0259 = IT_0123*IT_0215;
    const ccomplex_t IT_0260 = IT_0107*IT_0258*IT_0259;
    const ccomplex_t IT_0261 = 0.166666666666667*IT_0250 + -IT_0254 + -IT_0257
       + -IT_0260;
    const ccomplex_t IT_0262 = m_s*IT_0153;
    const ccomplex_t IT_0263 = (-0.5)*IT_0262;
    const ccomplex_t IT_0264 = 6*IT_0263;
    const ccomplex_t IT_0265 = IT_0009*IT_0143;
    const ccomplex_t IT_0266 = IT_0002*IT_0252*IT_0265;
    const ccomplex_t IT_0267 = IT_0109*IT_0225;
    const ccomplex_t IT_0268 = IT_0107*IT_0258*IT_0267;
    const ccomplex_t IT_0269 = m_s*IT_0146;
    const ccomplex_t IT_0270 = m_s*IT_0186;
    const ccomplex_t IT_0271 = m_s*IT_0207;
    const ccomplex_t IT_0272 = m_s*IT_0228;
    const ccomplex_t IT_0273 = 0.166666666666667*IT_0262 + IT_0266 + IT_0268 +
       -IT_0269 + -IT_0270 + -IT_0271 + -IT_0272;
    const ccomplex_t IT_0274 = 18*IT_0273;
    const ccomplex_t IT_0275 = IT_0040*IT_0183;
    const ccomplex_t IT_0276 = IT_0038*IT_0255*IT_0275;
    const ccomplex_t IT_0277 = IT_0073*IT_0204;
    const ccomplex_t IT_0278 = IT_0071*IT_0239*IT_0277;
    const ccomplex_t IT_0279 = -IT_0276 + -IT_0278;
    const ccomplex_t IT_0280 = (-18)*IT_0279;
    const ccomplex_t IT_0281 = IT_0264 + IT_0274 + IT_0280;
    const ccomplex_t IT_0282 = 6*conj(IT_0263);
    const ccomplex_t IT_0283 = 18*conj(IT_0273);
    const ccomplex_t IT_0284 = (-18)*conj(IT_0279);
    const ccomplex_t IT_0285 = IT_0282 + IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = conj(IT_0229)*(18*IT_0245 + (-6)*IT_0251 + (-18
      )*IT_0261) + 6*conj(IT_0230)*(IT_0245 + (-3)*IT_0251 + -IT_0261) + IT_0229
      *(18*conj(IT_0245) + (-6)*conj(IT_0251) + (-18)*conj(IT_0261)) + 6*IT_0230
      *(conj(IT_0245) + (-3)*conj(IT_0251) + -conj(IT_0261)) + conj(IT_0127)
      *IT_0281 + IT_0127*IT_0285;
    const ccomplex_t IT_0287 = m_b*m_s;
    const ccomplex_t IT_0288 = -conj(IT_0279);
    const ccomplex_t IT_0289 = 3*conj(IT_0263) + conj(IT_0273) + IT_0288;
    const ccomplex_t IT_0290 = 18*IT_0279;
    const ccomplex_t IT_0291 = (-18)*IT_0273;
    const ccomplex_t IT_0292 = (-6)*IT_0263;
    const ccomplex_t IT_0293 = IT_0290 + IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = 18*conj(IT_0279);
    const ccomplex_t IT_0295 = (-18)*conj(IT_0273);
    const ccomplex_t IT_0296 = (-6)*conj(IT_0263);
    const ccomplex_t IT_0297 = IT_0294 + IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = 18*IT_0287*(conj(IT_0251)*(IT_0263 + 1./3
      *IT_0273 + -1./3*IT_0279) + 0.0555555555555555*conj(IT_0261)*IT_0281 +
       0.0555555555555555*IT_0261*IT_0285 + 1./3*IT_0251*IT_0289 +
       0.0555555555555555*conj(IT_0245)*IT_0293 + 0.0555555555555555*IT_0245
      *IT_0297);
    const ccomplex_t IT_0299 = m_s*s_14;
    const ccomplex_t IT_0300 = conj(IT_0127)*((-18)*IT_0245 + 6*IT_0251 + 18
      *IT_0261) + IT_0127*((-18)*conj(IT_0245) + 6*conj(IT_0251) + 18*conj
      (IT_0261)) + (-18)*conj(IT_0230)*(IT_0263 + 1./3*IT_0273 + -1./3*IT_0279) 
      + (-6)*IT_0230*IT_0289 + conj(IT_0229)*IT_0293 + IT_0229*IT_0297;
    const ccomplex_t IT_0301 = 6*s_34;
    const ccomplex_t IT_0302 = 18*s_34;
    const ccomplex_t IT_0303 = (-18)*s_34;
    const ccomplex_t IT_0304 = (-6)*s_34;
    const ccomplex_t IT_0305 = 12*IT_0234;
    const ccomplex_t IT_0306 = (-6)*IT_0232;
    const ccomplex_t IT_0307 = (conj(IT_0127)*((-18)*IT_0229 + (-6)*IT_0230) +
       IT_0127*((-18)*conj(IT_0229) + (-6)*conj(IT_0230)))*IT_0231 + (IT_0127
      *conj(IT_0127) + IT_0229*conj(IT_0229) + IT_0230*conj(IT_0230))*IT_0236 +
       IT_0237*IT_0286 + IT_0298 + IT_0299*IT_0300 + (conj(IT_0251)*IT_0261 +
       IT_0251*conj(IT_0261) + conj(IT_0263)*IT_0273 + IT_0263*conj(IT_0273))
      *IT_0301 + (IT_0245*conj(IT_0245) + IT_0251*conj(IT_0251) + IT_0261*conj
      (IT_0261) + IT_0263*conj(IT_0263) + IT_0273*conj(IT_0273) + IT_0279*conj
      (IT_0279))*IT_0302 + (conj(IT_0245)*IT_0261 + IT_0245*conj(IT_0261) + conj
      (IT_0273)*IT_0279 + IT_0273*conj(IT_0279))*IT_0303 + (conj(IT_0245)
      *IT_0251 + IT_0245*conj(IT_0251) + conj(IT_0263)*IT_0279 + IT_0263*conj
      (IT_0279))*IT_0304 + (conj(IT_0229)*IT_0230 + IT_0229*conj(IT_0230))*
      (IT_0305 + IT_0306);
    return create_ccomplex_return(IT_0307);
}

