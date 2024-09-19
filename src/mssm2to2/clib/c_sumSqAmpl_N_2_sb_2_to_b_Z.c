#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_sb_2_to_b_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_sb_2_to_b_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_N_2, 2);
    const ccomplex_t IT_0002 = pow(m_sb_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*e_em*IT_0005*IT_0006;
    const ccomplex_t IT_0008 = 0.333333333333333*IT_0007;
    const ccomplex_t IT_0009 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0010 = IT_0005*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = cpow(IT_0006, -1);
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = m_b*N_d2*e_em*IT_0016*U_sb_01;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0011 + 1.5*IT_0019);
    const ccomplex_t IT_0021 = (-0.333333333333333)*IT_0020;
    const ccomplex_t IT_0022 = IT_0008*IT_0021;
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0005*IT_0006;
    const ccomplex_t IT_0027 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = IT_0004*IT_0014;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = IT_0029*IT_0031;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0028 + IT_0030 + 
      -IT_0032 + -IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0021*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0037*IT_0040;
    const ccomplex_t IT_0042 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0043 = IT_0026*IT_0042;
    const ccomplex_t IT_0044 = IT_0029*IT_0042;
    const ccomplex_t IT_0045 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0046 = IT_0026*IT_0045;
    const ccomplex_t IT_0047 = IT_0029*IT_0045;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0043 + IT_0044 + 
      -IT_0046 + -IT_0047);
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0051 = IT_0005*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = m_b*N_d1*e_em*IT_0016*U_sb_01;
    const ccomplex_t IT_0054 = IT_0015*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0052 + 1.5*IT_0055);
    const ccomplex_t IT_0057 = (-0.333333333333333)*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0049*IT_0060;
    const ccomplex_t IT_0062 = (-0.5)*IT_0041 + -IT_0061;
    const ccomplex_t IT_0063 = (-0.5)*IT_0034;
    const ccomplex_t IT_0064 = m_b*conj(N_d2)*e_em*IT_0016*U_sb_11;
    const ccomplex_t IT_0065 = IT_0015*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0068 = IT_0014*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0071 = IT_0005*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0066 + -IT_0069 +
       0.333333333333333*IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = (-3)*IT_0076;
    const ccomplex_t IT_0078 = 0.166666666666667*IT_0077;
    const ccomplex_t IT_0079 = IT_0038*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0063*IT_0080;
    const ccomplex_t IT_0082 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0083 = IT_0026*IT_0082;
    const ccomplex_t IT_0084 = IT_0029*IT_0082;
    const ccomplex_t IT_0085 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0086 = IT_0026*IT_0085;
    const ccomplex_t IT_0087 = IT_0029*IT_0085;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0083 + IT_0084 + 
      -IT_0086 + -IT_0087);
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = m_b*conj(N_d1)*e_em*IT_0016*U_sb_11;
    const ccomplex_t IT_0091 = IT_0015*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0094 = IT_0014*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0097 = IT_0005*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0092 + -IT_0095 +
       0.333333333333333*IT_0098);
    const ccomplex_t IT_0100 = 3*IT_0099;
    const ccomplex_t IT_0101 = 0.166666666666667*IT_0100;
    const ccomplex_t IT_0102 = IT_0058*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0089*IT_0103;
    const ccomplex_t IT_0105 = 0.5*IT_0081 + IT_0104;
    const ccomplex_t IT_0106 = e_em*IT_0026;
    const ccomplex_t IT_0107 = e_em*IT_0029;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0106 + 3*IT_0107);
    const ccomplex_t IT_0109 = (-0.166666666666667)*IT_0108;
    const ccomplex_t IT_0110 = IT_0078*IT_0109;
    const ccomplex_t IT_0111 = IT_0003*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = pow(m_Z, -2);
    const ccomplex_t IT_0114 = pow(m_Z, 4);
    const ccomplex_t IT_0115 = s_13*IT_0114;
    const ccomplex_t IT_0116 = IT_0113*IT_0115;
    const ccomplex_t IT_0117 = (-6)*IT_0116;
    const ccomplex_t IT_0118 = s_14*s_34;
    const ccomplex_t IT_0119 = 24*IT_0118;
    const ccomplex_t IT_0120 = IT_0117 + IT_0119;
    const ccomplex_t IT_0121 = pow(m_Z, 2);
    const ccomplex_t IT_0122 = s_13*IT_0121;
    const ccomplex_t IT_0123 = (-2)*IT_0118;
    const ccomplex_t IT_0124 = IT_0122 + IT_0123;
    const ccomplex_t IT_0125 = IT_0113*IT_0124;
    const ccomplex_t IT_0126 = (-6)*IT_0125;
    const ccomplex_t IT_0127 = 12*s_13;
    const ccomplex_t IT_0128 = IT_0126 + IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0130 = IT_0021*IT_0109;
    const ccomplex_t IT_0131 = IT_0003*IT_0129*IT_0130;
    const ccomplex_t IT_0132 = m_N_2*IT_0063;
    const ccomplex_t IT_0133 = 0.5*IT_0034;
    const ccomplex_t IT_0134 = m_N_2*IT_0133;
    const ccomplex_t IT_0135 = -IT_0134;
    const ccomplex_t IT_0136 = IT_0132 + IT_0135;
    const ccomplex_t IT_0137 = IT_0080*IT_0136;
    const ccomplex_t IT_0138 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0139 = IT_0014*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0142 = IT_0005*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = m_b*conj(N_d3)*e_em*IT_0016*U_sb_11;
    const ccomplex_t IT_0145 = IT_0015*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0140 + (
      -0.333333333333333)*IT_0143 + -IT_0146);
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = 3*IT_0148;
    const ccomplex_t IT_0150 = 0.166666666666667*IT_0149;
    const ccomplex_t IT_0151 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0152 = IT_0150*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0155 = IT_0026*IT_0154;
    const ccomplex_t IT_0156 = IT_0029*IT_0154;
    const ccomplex_t IT_0157 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0158 = IT_0026*IT_0157;
    const ccomplex_t IT_0159 = IT_0029*IT_0157;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*(IT_0155 + IT_0156 + 
      -IT_0158 + -IT_0159);
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = m_N_3*IT_0161;
    const ccomplex_t IT_0163 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0164 = IT_0026*IT_0163;
    const ccomplex_t IT_0165 = IT_0029*IT_0163;
    const ccomplex_t IT_0166 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0167 = IT_0026*IT_0166;
    const ccomplex_t IT_0168 = IT_0029*IT_0166;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(IT_0164 + IT_0165 + 
      -IT_0167 + -IT_0168);
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = m_N_2*IT_0170;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = IT_0162 + IT_0172;
    const ccomplex_t IT_0174 = IT_0153*IT_0173;
    const ccomplex_t IT_0175 = m_b*conj(N_d4)*e_em*IT_0016*U_sb_11;
    const ccomplex_t IT_0176 = IT_0015*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0179 = IT_0005*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0182 = IT_0014*IT_0181;
    const ccomplex_t IT_0183 = 1.4142135623731*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*(IT_0177 + 0.333333333333333
      *IT_0180 + -IT_0183);
    const ccomplex_t IT_0185 = 3*IT_0184;
    const ccomplex_t IT_0186 = 0.166666666666667*IT_0185;
    const ccomplex_t IT_0187 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0188 = IT_0186*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0191 = IT_0026*IT_0190;
    const ccomplex_t IT_0192 = IT_0029*IT_0190;
    const ccomplex_t IT_0193 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0194 = IT_0026*IT_0193;
    const ccomplex_t IT_0195 = IT_0029*IT_0193;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0191 + IT_0192 + 
      -IT_0194 + -IT_0195);
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = m_N_4*IT_0197;
    const ccomplex_t IT_0199 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0200 = IT_0029*IT_0199;
    const ccomplex_t IT_0201 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0202 = IT_0026*IT_0201;
    const ccomplex_t IT_0203 = IT_0029*IT_0201;
    const ccomplex_t IT_0204 = IT_0026*IT_0199;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*(IT_0200 + -IT_0202 + 
      -IT_0203 + IT_0204);
    const ccomplex_t IT_0206 = -IT_0205;
    const ccomplex_t IT_0207 = (-0.5)*IT_0206;
    const ccomplex_t IT_0208 = m_N_2*IT_0207;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = IT_0198 + IT_0209;
    const ccomplex_t IT_0211 = IT_0189*IT_0210;
    const ccomplex_t IT_0212 = m_b*IT_0024;
    const ccomplex_t IT_0213 = -IT_0131 + (-0.5)*IT_0137 + -IT_0174 + -IT_0211
       + IT_0212;
    const ccomplex_t IT_0214 = m_N_1*IT_0049;
    const ccomplex_t IT_0215 = m_N_2*IT_0089;
    const ccomplex_t IT_0216 = -IT_0215;
    const ccomplex_t IT_0217 = IT_0214 + IT_0216;
    const ccomplex_t IT_0218 = IT_0103*IT_0217;
    const ccomplex_t IT_0219 = (-0.5)*IT_0137 + IT_0218;
    const ccomplex_t IT_0220 = IT_0213 + IT_0219;
    const ccomplex_t IT_0221 = conj(IT_0213) + conj(IT_0219);
    const ccomplex_t IT_0222 = IT_0008*IT_0078;
    const ccomplex_t IT_0223 = IT_0003*IT_0129*IT_0222;
    const ccomplex_t IT_0224 = m_b*IT_0112;
    const ccomplex_t IT_0225 = IT_0040*IT_0136;
    const ccomplex_t IT_0226 = m_b*N_d3*e_em*IT_0016*U_sb_01;
    const ccomplex_t IT_0227 = IT_0015*IT_0226;
    const ccomplex_t IT_0228 = 1.4142135623731*IT_0227;
    const ccomplex_t IT_0229 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0230 = IT_0005*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*(IT_0228 + 0.666666666666667
      *IT_0231);
    const ccomplex_t IT_0233 = 1.5*IT_0232;
    const ccomplex_t IT_0234 = (-0.333333333333333)*IT_0233;
    const ccomplex_t IT_0235 = IT_0151*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = m_N_3*IT_0170;
    const ccomplex_t IT_0238 = m_N_2*IT_0161;
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = IT_0237 + IT_0239;
    const ccomplex_t IT_0241 = IT_0236*IT_0240;
    const ccomplex_t IT_0242 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0243 = IT_0005*IT_0242;
    const ccomplex_t IT_0244 = 1.4142135623731*IT_0243;
    const ccomplex_t IT_0245 = m_b*N_d4*e_em*IT_0016*U_sb_01;
    const ccomplex_t IT_0246 = IT_0015*IT_0245;
    const ccomplex_t IT_0247 = 1.4142135623731*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*(IT_0244 + 1.5*IT_0247);
    const ccomplex_t IT_0249 = (-0.333333333333333)*IT_0248;
    const ccomplex_t IT_0250 = IT_0187*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*-1)*(IT_0200 + -IT_0202 + 
      -IT_0203 + IT_0204);
    const ccomplex_t IT_0253 = (-0.5)*IT_0252;
    const ccomplex_t IT_0254 = m_N_4*IT_0253;
    const ccomplex_t IT_0255 = m_N_2*IT_0197;
    const ccomplex_t IT_0256 = -IT_0255;
    const ccomplex_t IT_0257 = IT_0254 + IT_0256;
    const ccomplex_t IT_0258 = IT_0251*IT_0257;
    const ccomplex_t IT_0259 = IT_0223 + -IT_0224 + (-0.5)*IT_0225 + IT_0241 +
       IT_0258;
    const ccomplex_t IT_0260 = m_N_1*IT_0089;
    const ccomplex_t IT_0261 = m_N_2*IT_0049;
    const ccomplex_t IT_0262 = -IT_0261;
    const ccomplex_t IT_0263 = IT_0260 + IT_0262;
    const ccomplex_t IT_0264 = IT_0060*IT_0263;
    const ccomplex_t IT_0265 = (-0.5)*IT_0225 + -IT_0264;
    const ccomplex_t IT_0266 = IT_0259 + IT_0265;
    const ccomplex_t IT_0267 = conj(IT_0259) + conj(IT_0265);
    const ccomplex_t IT_0268 = (-24)*IT_0118;
    const ccomplex_t IT_0269 = 6*IT_0116;
    const ccomplex_t IT_0270 = IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = IT_0037*IT_0080;
    const ccomplex_t IT_0272 = (0 + _Complex_I*-1)*(IT_0140 + (
      -0.333333333333333)*IT_0143 + -IT_0146);
    const ccomplex_t IT_0273 = 3*IT_0272;
    const ccomplex_t IT_0274 = 0.166666666666667*IT_0273;
    const ccomplex_t IT_0275 = IT_0151*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*IT_0275;
    const ccomplex_t IT_0277 = IT_0170*IT_0276;
    const ccomplex_t IT_0278 = IT_0189*IT_0207;
    const ccomplex_t IT_0279 = 0.5*IT_0271 + IT_0277 + IT_0278;
    const ccomplex_t IT_0280 = IT_0161*IT_0236;
    const ccomplex_t IT_0281 = IT_0197*IT_0251;
    const ccomplex_t IT_0282 = -IT_0280 + -IT_0281;
    const ccomplex_t IT_0283 = IT_0040*IT_0063;
    const ccomplex_t IT_0284 = (-0.5)*IT_0283;
    const ccomplex_t IT_0285 = IT_0282 + IT_0284;
    const ccomplex_t IT_0286 = conj(IT_0282) + conj(IT_0284);
    const ccomplex_t IT_0287 = (-12)*IT_0122;
    const ccomplex_t IT_0288 = IT_0119 + IT_0269 + IT_0287;
    const ccomplex_t IT_0289 = 12*IT_0122;
    const ccomplex_t IT_0290 = IT_0117 + IT_0268 + IT_0289;
    const ccomplex_t IT_0291 = m_b*m_N_2;
    const ccomplex_t IT_0292 = pow(s_14, 2);
    const ccomplex_t IT_0293 = IT_0113*IT_0292;
    const ccomplex_t IT_0294 = -IT_0001;
    const ccomplex_t IT_0295 = IT_0293 + IT_0294;
    const ccomplex_t IT_0296 = IT_0291*IT_0295;
    const ccomplex_t IT_0297 = IT_0080*IT_0133;
    const ccomplex_t IT_0298 = IT_0153*IT_0170;
    const ccomplex_t IT_0299 = IT_0189*IT_0253;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0005
      *IT_0006 + 3*IT_0004*IT_0014)*U_sb_01*conj(U_sb_01) + IT_0005*IT_0006
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0301 = (-0.666666666666667)*IT_0300;
    const ccomplex_t IT_0302 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0303 = IT_0078*IT_0302;
    const ccomplex_t IT_0304 = (0 + _Complex_I*1)*IT_0303;
    const ccomplex_t IT_0305 = IT_0301*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0005
      *IT_0006 + 3*IT_0004*IT_0014)*conj(U_sb_00)*U_sb_01 + IT_0005*IT_0006*conj
      (U_sb_10)*U_sb_11);
    const ccomplex_t IT_0307 = (-0.666666666666667)*IT_0306;
    const ccomplex_t IT_0308 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0309 = IT_0014*IT_0308;
    const ccomplex_t IT_0310 = 1.4142135623731*IT_0309;
    const ccomplex_t IT_0311 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0312 = IT_0005*IT_0311;
    const ccomplex_t IT_0313 = 1.4142135623731*IT_0312;
    const ccomplex_t IT_0314 = m_b*conj(N_d2)*e_em*IT_0016*U_sb_10;
    const ccomplex_t IT_0315 = IT_0015*IT_0314;
    const ccomplex_t IT_0316 = 1.4142135623731*IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*(IT_0310 + (
      -0.333333333333333)*IT_0313 + -IT_0316);
    const ccomplex_t IT_0318 = -IT_0317;
    const ccomplex_t IT_0319 = -IT_0318;
    const ccomplex_t IT_0320 = -IT_0319;
    const ccomplex_t IT_0321 = 3*IT_0320;
    const ccomplex_t IT_0322 = 0.333333333333333*IT_0321;
    const ccomplex_t IT_0323 = 3*IT_0322;
    const ccomplex_t IT_0324 = (-0.333333333333333)*IT_0323;
    const ccomplex_t IT_0325 = -IT_0324;
    const ccomplex_t IT_0326 = 3*IT_0325;
    const ccomplex_t IT_0327 = 0.166666666666667*IT_0326;
    const ccomplex_t IT_0328 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0329 = IT_0327*IT_0328;
    const ccomplex_t IT_0330 = (0 + _Complex_I*1)*IT_0329;
    const ccomplex_t IT_0331 = IT_0307*IT_0330;
    const ccomplex_t IT_0332 = -IT_0297 + (-2)*IT_0298 + (-2)*IT_0299 + 
      -IT_0305 + -IT_0331;
    const ccomplex_t IT_0333 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0334 = IT_0005*IT_0333;
    const ccomplex_t IT_0335 = 1.4142135623731*IT_0334;
    const ccomplex_t IT_0336 = m_b*N_d2*e_em*IT_0016*U_sb_00;
    const ccomplex_t IT_0337 = IT_0015*IT_0336;
    const ccomplex_t IT_0338 = 1.4142135623731*IT_0337;
    const ccomplex_t IT_0339 = (0 + _Complex_I*1)*(IT_0335 + 1.5*IT_0338);
    const ccomplex_t IT_0340 = 0.666666666666667*IT_0339;
    const ccomplex_t IT_0341 = 1.5*IT_0340;
    const ccomplex_t IT_0342 = 0.666666666666667*IT_0341;
    const ccomplex_t IT_0343 = 1.5*IT_0342;
    const ccomplex_t IT_0344 = 0.666666666666667*IT_0343;
    const ccomplex_t IT_0345 = 1.5*IT_0344;
    const ccomplex_t IT_0346 = 0.666666666666667*IT_0345;
    const ccomplex_t IT_0347 = 1.5*IT_0346;
    const ccomplex_t IT_0348 = (-0.333333333333333)*IT_0347;
    const ccomplex_t IT_0349 = IT_0328*IT_0348;
    const ccomplex_t IT_0350 = (0 + _Complex_I*1)*IT_0349;
    const ccomplex_t IT_0351 = IT_0307*IT_0350;
    const ccomplex_t IT_0352 = (0 + _Complex_I*1)*(1.5*IT_0228 + IT_0231);
    const ccomplex_t IT_0353 = (-0.333333333333333)*IT_0352;
    const ccomplex_t IT_0354 = IT_0151*IT_0353;
    const ccomplex_t IT_0355 = (0 + _Complex_I*1)*IT_0354;
    const ccomplex_t IT_0356 = IT_0161*IT_0355;
    const ccomplex_t IT_0357 = IT_0021*IT_0302;
    const ccomplex_t IT_0358 = (0 + _Complex_I*1)*IT_0357;
    const ccomplex_t IT_0359 = IT_0301*IT_0358;
    const ccomplex_t IT_0360 = 2*IT_0281 + IT_0283 + IT_0351 + 2*IT_0356 +
       IT_0359;
    const ccomplex_t IT_0361 = 6*IT_0360;
    const ccomplex_t IT_0362 = IT_0040*IT_0133;
    const ccomplex_t IT_0363 = 2*IT_0061 + IT_0362;
    const ccomplex_t IT_0364 = (-6)*IT_0363;
    const ccomplex_t IT_0365 = IT_0361 + IT_0364;
    const ccomplex_t IT_0366 = 6*conj(IT_0360);
    const ccomplex_t IT_0367 = (-6)*conj(IT_0363);
    const ccomplex_t IT_0368 = IT_0366 + IT_0367;
    const ccomplex_t IT_0369 = -IT_0081 + (-2)*IT_0104;
    const ccomplex_t IT_0370 = 6*IT_0363;
    const ccomplex_t IT_0371 = (-6)*IT_0360;
    const ccomplex_t IT_0372 = IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = 6*conj(IT_0363);
    const ccomplex_t IT_0374 = (-6)*conj(IT_0360);
    const ccomplex_t IT_0375 = IT_0373 + IT_0374;
    const ccomplex_t IT_0376 = m_b*IT_0295;
    const ccomplex_t IT_0377 = 6*IT_0332;
    const ccomplex_t IT_0378 = (-6)*IT_0369;
    const ccomplex_t IT_0379 = IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = 6*conj(IT_0332);
    const ccomplex_t IT_0381 = (-6)*conj(IT_0369);
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = s_14*s_34*IT_0113;
    const ccomplex_t IT_0384 = -IT_0383;
    const ccomplex_t IT_0385 = s_13 + IT_0384;
    const ccomplex_t IT_0386 = s_13*IT_0385;
    const ccomplex_t IT_0387 = (-2)*IT_0024 + -IT_0351 + -IT_0359;
    const ccomplex_t IT_0388 = 2*IT_0112 + IT_0305 + IT_0331;
    const ccomplex_t IT_0389 = 6*IT_0369;
    const ccomplex_t IT_0390 = (-6)*IT_0332;
    const ccomplex_t IT_0391 = IT_0389 + IT_0390;
    const ccomplex_t IT_0392 = 6*conj(IT_0369);
    const ccomplex_t IT_0393 = (-6)*conj(IT_0332);
    const ccomplex_t IT_0394 = IT_0392 + IT_0393;
    const ccomplex_t IT_0395 = s_13*s_14;
    const ccomplex_t IT_0396 = s_34*IT_0001;
    const ccomplex_t IT_0397 = -IT_0396;
    const ccomplex_t IT_0398 = IT_0395 + IT_0397;
    const ccomplex_t IT_0399 = conj(IT_0025) + conj(IT_0062);
    const ccomplex_t IT_0400 = IT_0025 + IT_0062;
    const ccomplex_t IT_0401 = conj(IT_0105) + conj(IT_0112);
    const ccomplex_t IT_0402 = IT_0105 + IT_0112;
    const ccomplex_t IT_0403 = s_34*IT_0121;
    const ccomplex_t IT_0404 = IT_0113*IT_0403;
    const ccomplex_t IT_0405 = 0.5*IT_0404;
    const ccomplex_t IT_0406 = s_34 + IT_0405;
    const ccomplex_t IT_0407 = m_N_2*IT_0406;
    const ccomplex_t IT_0408 = 2*IT_0407;
    const ccomplex_t IT_0409 = 6*IT_0213;
    const ccomplex_t IT_0410 = 6*IT_0219;
    const ccomplex_t IT_0411 = IT_0409 + IT_0410;
    const ccomplex_t IT_0412 = 6*conj(IT_0213);
    const ccomplex_t IT_0413 = 6*conj(IT_0219);
    const ccomplex_t IT_0414 = IT_0412 + IT_0413;
    const ccomplex_t IT_0415 = 6*IT_0259;
    const ccomplex_t IT_0416 = 6*IT_0265;
    const ccomplex_t IT_0417 = IT_0415 + IT_0416;
    const ccomplex_t IT_0418 = 6*conj(IT_0259);
    const ccomplex_t IT_0419 = 6*conj(IT_0265);
    const ccomplex_t IT_0420 = IT_0418 + IT_0419;
    const ccomplex_t IT_0421 = s_14*IT_0121;
    const ccomplex_t IT_0422 = IT_0113*IT_0421;
    const ccomplex_t IT_0423 = (-0.25)*IT_0422;
    const ccomplex_t IT_0424 = s_14 + IT_0423;
    const ccomplex_t IT_0425 = m_b*IT_0424;
    const ccomplex_t IT_0426 = (-4)*IT_0425;
    const ccomplex_t IT_0427 = (-0.25)*IT_0404;
    const ccomplex_t IT_0428 = s_34 + IT_0427;
    const ccomplex_t IT_0429 = m_N_2*IT_0428;
    const ccomplex_t IT_0430 = (-4)*IT_0429;
    const ccomplex_t IT_0431 = (-6)*IT_0213;
    const ccomplex_t IT_0432 = (-6)*IT_0219;
    const ccomplex_t IT_0433 = IT_0431 + IT_0432;
    const ccomplex_t IT_0434 = (-6)*conj(IT_0213);
    const ccomplex_t IT_0435 = (-6)*conj(IT_0219);
    const ccomplex_t IT_0436 = IT_0434 + IT_0435;
    const ccomplex_t IT_0437 = (-6)*IT_0259;
    const ccomplex_t IT_0438 = (-6)*IT_0265;
    const ccomplex_t IT_0439 = IT_0437 + IT_0438;
    const ccomplex_t IT_0440 = (-6)*conj(IT_0259);
    const ccomplex_t IT_0441 = (-6)*conj(IT_0265);
    const ccomplex_t IT_0442 = IT_0440 + IT_0441;
    const ccomplex_t IT_0443 = 0.5*IT_0422;
    const ccomplex_t IT_0444 = s_14 + IT_0443;
    const ccomplex_t IT_0445 = m_b*IT_0444;
    const ccomplex_t IT_0446 = 2*IT_0445;
    const ccomplex_t IT_0447 = IT_0121*IT_0291;
    const ccomplex_t IT_0448 = 3*IT_0447;
    const ccomplex_t IT_0449 = 6*IT_0284;
    const ccomplex_t IT_0450 = (-0.166666666666667)*IT_0449;
    const ccomplex_t IT_0451 = 6*IT_0282;
    const ccomplex_t IT_0452 = (-0.166666666666667)*IT_0451;
    const ccomplex_t IT_0453 = IT_0062 + IT_0450 + IT_0452;
    const ccomplex_t IT_0454 = 6*conj(IT_0284);
    const ccomplex_t IT_0455 = (-0.166666666666667)*IT_0454;
    const ccomplex_t IT_0456 = 6*conj(IT_0282);
    const ccomplex_t IT_0457 = (-0.166666666666667)*IT_0456;
    const ccomplex_t IT_0458 = conj(IT_0062) + IT_0455 + IT_0457;
    const ccomplex_t IT_0459 = 6*IT_0025;
    const ccomplex_t IT_0460 = 6*conj(IT_0025);
    const ccomplex_t IT_0461 = (-3)*IT_0447;
    const ccomplex_t IT_0462 = 6*(conj(IT_0105)*(IT_0062 + -IT_0285) + IT_0105
      *(conj(IT_0062) + -IT_0286) + -conj(IT_0279)*IT_0453 + -IT_0279*IT_0458 +
       0.166666666666667*conj(IT_0112)*IT_0459 + 0.166666666666667*IT_0112
      *IT_0460)*IT_0461;
    const ccomplex_t IT_0463 = (-18)*conj(IT_0213);
    const ccomplex_t IT_0464 = s_13*IT_0295;
    const ccomplex_t IT_0465 = IT_0291*IT_0385;
    const ccomplex_t IT_0466 = 6*IT_0387;
    const ccomplex_t IT_0467 = 6*conj(IT_0387);
    const ccomplex_t IT_0468 = (-6)*IT_0387;
    const ccomplex_t IT_0469 = (-6)*conj(IT_0387);
    const ccomplex_t IT_0470 = pow(s_34, 2);
    const ccomplex_t IT_0471 = IT_0113*IT_0470;
    const ccomplex_t IT_0472 = -IT_0471;
    const ccomplex_t IT_0473 = IT_0000 + IT_0472;
    const ccomplex_t IT_0474 = IT_0291*IT_0473;
    const ccomplex_t IT_0475 = s_13*IT_0473;
    const ccomplex_t IT_0476 = (-6)*IT_0388;
    const ccomplex_t IT_0477 = m_b*IT_0385;
    const ccomplex_t IT_0478 = (-6)*conj(IT_0388);
    const ccomplex_t IT_0479 = s_14*IT_0000;
    const ccomplex_t IT_0480 = s_13*s_34;
    const ccomplex_t IT_0481 = -IT_0480;
    const ccomplex_t IT_0482 = IT_0479 + IT_0481;
    const ccomplex_t IT_0483 = conj(IT_0279)*IT_0388 + IT_0279*conj(IT_0388) +
       0.166666666666667*IT_0286*IT_0466 + 0.166666666666667*IT_0285*IT_0467 +
       0.166666666666667*IT_0399*IT_0468 + 0.166666666666667*IT_0400*IT_0469 +
       0.166666666666667*IT_0401*IT_0476 + 0.166666666666667*IT_0402*IT_0478;
    const ccomplex_t IT_0484 = m_N_2*IT_0385;
    const ccomplex_t IT_0485 = m_N_2*IT_0473;
    const ccomplex_t IT_0486 = (conj(IT_0025)*IT_0062 + IT_0025*conj(IT_0062) 
      + conj(IT_0105)*IT_0112 + IT_0105*conj(IT_0112))*IT_0120 + IT_0128*
      (IT_0220*IT_0221 + IT_0266*IT_0267) + IT_0270*(conj(IT_0112)*IT_0279 +
       IT_0112*conj(IT_0279) + conj(IT_0025)*IT_0285 + IT_0025*IT_0286) + 
      (IT_0025*conj(IT_0025) + IT_0062*conj(IT_0062) + IT_0105*conj(IT_0105) +
       IT_0112*conj(IT_0112) + IT_0279*conj(IT_0279) + (IT_0282 + IT_0284)
      *IT_0286)*IT_0288 + (conj(IT_0105)*IT_0279 + IT_0105*conj(IT_0279) + conj
      (IT_0062)*IT_0285 + IT_0062*IT_0286)*IT_0290 + IT_0296*(conj(IT_0332)
      *IT_0365 + IT_0332*IT_0368 + conj(IT_0369)*IT_0372 + IT_0369*IT_0375) +
       IT_0376*(IT_0221*IT_0365 + IT_0220*IT_0368 + IT_0267*IT_0379 + IT_0266
      *IT_0382) + IT_0386*(IT_0375*IT_0387 + IT_0372*conj(IT_0387) + conj
      (IT_0388)*IT_0391 + IT_0388*IT_0394) + IT_0398*(IT_0286*IT_0365 + IT_0285
      *IT_0368 + conj(IT_0279)*IT_0379 + IT_0279*IT_0382 + IT_0372*IT_0399 +
       IT_0375*IT_0400 + IT_0391*IT_0401 + IT_0394*IT_0402) + IT_0408*(conj
      (IT_0112)*IT_0411 + IT_0112*IT_0414 + conj(IT_0025)*IT_0417 + IT_0025
      *IT_0420) + (conj(IT_0025)*IT_0411 + IT_0025*IT_0414 + conj(IT_0112)
      *IT_0417 + IT_0112*IT_0420)*IT_0426 + IT_0430*(conj(IT_0279)*IT_0411 +
       IT_0279*IT_0414 + IT_0286*IT_0417 + IT_0285*IT_0420 + conj(IT_0105)
      *IT_0433 + IT_0105*IT_0436 + conj(IT_0062)*IT_0439 + IT_0062*IT_0442) + 
      (IT_0286*IT_0411 + IT_0285*IT_0414 + conj(IT_0279)*IT_0417 + IT_0279
      *IT_0420 + conj(IT_0062)*IT_0433 + IT_0062*IT_0436 + conj(IT_0105)*IT_0439
       + IT_0105*IT_0442)*IT_0446 + (-6)*IT_0448*(conj(IT_0025)*IT_0105 +
       IT_0025*conj(IT_0105) + conj(IT_0112)*IT_0453 + IT_0112*IT_0458 + (
      -0.166666666666667)*conj(IT_0279)*IT_0459 + (-0.166666666666667)*IT_0279
      *IT_0460) + IT_0462 + IT_0291*(IT_0220*((-18)*conj(IT_0259) + (-18)*conj
      (IT_0265)) + IT_0266*((-18)*conj(IT_0219) + IT_0463)) + (conj(IT_0360)
      *IT_0361 + conj(IT_0363)*IT_0372 + IT_0363*IT_0374 + conj(IT_0332)*IT_0377
       + conj(IT_0369)*IT_0391 + IT_0369*IT_0393)*IT_0464 + IT_0465*(IT_0375
      *IT_0388 + IT_0372*conj(IT_0388) + conj(IT_0369)*IT_0466 + IT_0369*IT_0467
       + conj(IT_0332)*IT_0468 + IT_0332*IT_0469) + (conj(IT_0388)*IT_0468 +
       IT_0388*IT_0469)*IT_0474 + IT_0475*(conj(IT_0387)*IT_0468 + conj(IT_0388)
      *IT_0476) + IT_0477*(IT_0221*IT_0468 + IT_0220*IT_0469 + IT_0267*IT_0476 +
       IT_0266*IT_0478) + 6*IT_0482*IT_0483 + (IT_0267*IT_0372 + IT_0266*IT_0375
       + conj(IT_0219)*IT_0391 + IT_0219*IT_0394 + conj(IT_0369)*IT_0409 +
       IT_0369*IT_0412 + conj(IT_0332)*IT_0431 + IT_0332*IT_0434)*IT_0484 + 
      (conj(IT_0388)*IT_0431 + IT_0388*IT_0434 + IT_0267*IT_0468 + IT_0266
      *IT_0469 + conj(IT_0219)*IT_0476 + IT_0219*IT_0478)*IT_0485;
    return create_ccomplex_return(IT_0486);
}

