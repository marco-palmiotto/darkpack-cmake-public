#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_sb_1_to_b_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_sb_1_to_b_Z(
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
    const ccomplex_t IT_0001 = pow(m_N_1, 2);
    const ccomplex_t IT_0002 = pow(m_sb_1, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*e_em*IT_0005*IT_0006;
    const ccomplex_t IT_0008 = 0.333333333333333*IT_0007;
    const ccomplex_t IT_0009 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0010 = IT_0005*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = cpow(IT_0006, -1);
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = m_b*N_d1*e_em*IT_0016*U_sb_00;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0011 + 1.5*IT_0019);
    const ccomplex_t IT_0021 = (-0.333333333333333)*IT_0020;
    const ccomplex_t IT_0022 = IT_0008*IT_0021;
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0005*IT_0006;
    const ccomplex_t IT_0027 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = IT_0004*IT_0014;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = IT_0029*IT_0031;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0028 + IT_0030 + 
      -IT_0032 + -IT_0033);
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0021*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0035*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = (-0.5)*IT_0034;
    const ccomplex_t IT_0042 = m_b*conj(N_d1)*e_em*IT_0016*U_sb_10;
    const ccomplex_t IT_0043 = IT_0015*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0046 = IT_0005*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0049 = IT_0014*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0044 + 0.333333333333333
      *IT_0047 + -IT_0050);
    const ccomplex_t IT_0052 = 3*IT_0051;
    const ccomplex_t IT_0053 = 0.166666666666667*IT_0052;
    const ccomplex_t IT_0054 = IT_0036*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0041*IT_0055;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = e_em*IT_0026;
    const ccomplex_t IT_0059 = e_em*IT_0029;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0058 + 3*IT_0059);
    const ccomplex_t IT_0061 = (-0.166666666666667)*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = pow(m_Z, -2);
    const ccomplex_t IT_0066 = pow(m_Z, 4);
    const ccomplex_t IT_0067 = s_13*IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = (-6)*IT_0068;
    const ccomplex_t IT_0070 = s_14*s_34;
    const ccomplex_t IT_0071 = 24*IT_0070;
    const ccomplex_t IT_0072 = IT_0069 + IT_0071;
    const ccomplex_t IT_0073 = pow(m_Z, 2);
    const ccomplex_t IT_0074 = s_13*IT_0073;
    const ccomplex_t IT_0075 = (-2)*IT_0070;
    const ccomplex_t IT_0076 = IT_0074 + IT_0075;
    const ccomplex_t IT_0077 = IT_0065*IT_0076;
    const ccomplex_t IT_0078 = (-6)*IT_0077;
    const ccomplex_t IT_0079 = 12*s_13;
    const ccomplex_t IT_0080 = IT_0078 + IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0082 = IT_0021*IT_0061;
    const ccomplex_t IT_0083 = IT_0003*IT_0081*IT_0082;
    const ccomplex_t IT_0084 = m_b*IT_0024;
    const ccomplex_t IT_0085 = m_N_1*IT_0041;
    const ccomplex_t IT_0086 = m_N_1*IT_0035;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0085 + IT_0087;
    const ccomplex_t IT_0089 = IT_0055*IT_0088;
    const ccomplex_t IT_0090 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0091 = IT_0026*IT_0090;
    const ccomplex_t IT_0092 = IT_0029*IT_0090;
    const ccomplex_t IT_0093 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0094 = IT_0026*IT_0093;
    const ccomplex_t IT_0095 = IT_0029*IT_0093;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0091 + IT_0092 + 
      -IT_0094 + -IT_0095);
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = m_N_2*IT_0097;
    const ccomplex_t IT_0099 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0100 = IT_0026*IT_0099;
    const ccomplex_t IT_0101 = IT_0029*IT_0099;
    const ccomplex_t IT_0102 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0103 = IT_0026*IT_0102;
    const ccomplex_t IT_0104 = IT_0029*IT_0102;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0100 + IT_0101 + 
      -IT_0103 + -IT_0104);
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = m_N_1*IT_0106;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = IT_0098 + IT_0108;
    const ccomplex_t IT_0110 = m_b*conj(N_d2)*e_em*IT_0016*U_sb_10;
    const ccomplex_t IT_0111 = IT_0015*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0114 = IT_0014*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0117 = IT_0005*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0112 + -IT_0115 +
       0.333333333333333*IT_0118);
    const ccomplex_t IT_0120 = 3*IT_0119;
    const ccomplex_t IT_0121 = 0.166666666666667*IT_0120;
    const ccomplex_t IT_0122 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0123 = IT_0121*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = IT_0109*IT_0124;
    const ccomplex_t IT_0126 = m_b*conj(N_d3)*e_em*IT_0016*U_sb_10;
    const ccomplex_t IT_0127 = IT_0015*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0130 = IT_0005*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0133 = IT_0014*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0128 + 0.333333333333333
      *IT_0131 + -IT_0134);
    const ccomplex_t IT_0136 = 3*IT_0135;
    const ccomplex_t IT_0137 = 0.166666666666667*IT_0136;
    const ccomplex_t IT_0138 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0142 = IT_0026*IT_0141;
    const ccomplex_t IT_0143 = IT_0029*IT_0141;
    const ccomplex_t IT_0144 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0145 = IT_0026*IT_0144;
    const ccomplex_t IT_0146 = IT_0029*IT_0144;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0142 + IT_0143 + 
      -IT_0145 + -IT_0146);
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = m_N_3*IT_0148;
    const ccomplex_t IT_0150 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0151 = IT_0026*IT_0150;
    const ccomplex_t IT_0152 = IT_0029*IT_0150;
    const ccomplex_t IT_0153 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0154 = IT_0026*IT_0153;
    const ccomplex_t IT_0155 = IT_0029*IT_0153;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0151 + IT_0152 + 
      -IT_0154 + -IT_0155);
    const ccomplex_t IT_0157 = (-0.5)*IT_0156;
    const ccomplex_t IT_0158 = m_N_1*IT_0157;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = IT_0149 + IT_0159;
    const ccomplex_t IT_0161 = IT_0140*IT_0160;
    const ccomplex_t IT_0162 = m_b*conj(N_d4)*e_em*IT_0016*U_sb_10;
    const ccomplex_t IT_0163 = IT_0015*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0166 = IT_0014*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0169 = IT_0005*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*(IT_0164 + -IT_0167 +
       0.333333333333333*IT_0170);
    const ccomplex_t IT_0172 = 3*IT_0171;
    const ccomplex_t IT_0173 = 0.166666666666667*IT_0172;
    const ccomplex_t IT_0174 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0175 = IT_0173*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0178 = IT_0026*IT_0177;
    const ccomplex_t IT_0179 = IT_0029*IT_0177;
    const ccomplex_t IT_0180 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0181 = IT_0026*IT_0180;
    const ccomplex_t IT_0182 = IT_0029*IT_0180;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0178 + IT_0179 + 
      -IT_0181 + -IT_0182);
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = m_N_4*IT_0184;
    const ccomplex_t IT_0186 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0187 = IT_0026*IT_0186;
    const ccomplex_t IT_0188 = IT_0029*IT_0186;
    const ccomplex_t IT_0189 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0190 = IT_0026*IT_0189;
    const ccomplex_t IT_0191 = IT_0029*IT_0189;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*(IT_0187 + IT_0188 + 
      -IT_0190 + -IT_0191);
    const ccomplex_t IT_0193 = (-0.5)*IT_0192;
    const ccomplex_t IT_0194 = m_N_1*IT_0193;
    const ccomplex_t IT_0195 = -IT_0194;
    const ccomplex_t IT_0196 = IT_0185 + IT_0195;
    const ccomplex_t IT_0197 = IT_0176*IT_0196;
    const ccomplex_t IT_0198 = -IT_0083 + IT_0084 + (-0.5)*IT_0089 + -IT_0125 
      + -IT_0161 + -IT_0197;
    const ccomplex_t IT_0199 = IT_0055*(IT_0085 + -IT_0086);
    const ccomplex_t IT_0200 = (-0.5)*IT_0199;
    const ccomplex_t IT_0201 = IT_0198 + IT_0200;
    const ccomplex_t IT_0202 = conj(IT_0198) + conj(IT_0200);
    const ccomplex_t IT_0203 = m_N_3*IT_0157;
    const ccomplex_t IT_0204 = m_N_1*IT_0148;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = IT_0203 + IT_0205;
    const ccomplex_t IT_0207 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0208 = IT_0005*IT_0207;
    const ccomplex_t IT_0209 = 1.4142135623731*IT_0208;
    const ccomplex_t IT_0210 = m_b*N_d3*e_em*IT_0016*U_sb_00;
    const ccomplex_t IT_0211 = IT_0015*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*(IT_0209 + 1.5*IT_0212);
    const ccomplex_t IT_0214 = (-0.333333333333333)*IT_0213;
    const ccomplex_t IT_0215 = IT_0138*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = IT_0206*IT_0216;
    const ccomplex_t IT_0218 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0219 = IT_0005*IT_0218;
    const ccomplex_t IT_0220 = 1.4142135623731*IT_0219;
    const ccomplex_t IT_0221 = m_b*N_d4*e_em*IT_0016*U_sb_00;
    const ccomplex_t IT_0222 = IT_0015*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*(IT_0220 + 1.5*IT_0223);
    const ccomplex_t IT_0225 = (-0.333333333333333)*IT_0224;
    const ccomplex_t IT_0226 = IT_0174*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = m_N_4*IT_0193;
    const ccomplex_t IT_0229 = m_N_1*IT_0184;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = IT_0227*IT_0231;
    const ccomplex_t IT_0233 = IT_0008*IT_0053;
    const ccomplex_t IT_0234 = IT_0003*IT_0081*IT_0233;
    const ccomplex_t IT_0235 = m_b*IT_0064;
    const ccomplex_t IT_0236 = IT_0038*IT_0088;
    const ccomplex_t IT_0237 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0238 = IT_0005*IT_0237;
    const ccomplex_t IT_0239 = 1.4142135623731*IT_0238;
    const ccomplex_t IT_0240 = m_b*N_d2*e_em*IT_0016*U_sb_00;
    const ccomplex_t IT_0241 = IT_0015*IT_0240;
    const ccomplex_t IT_0242 = 1.4142135623731*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*(IT_0239 + 1.5*IT_0242);
    const ccomplex_t IT_0244 = (-0.333333333333333)*IT_0243;
    const ccomplex_t IT_0245 = IT_0122*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*IT_0245;
    const ccomplex_t IT_0247 = m_N_2*IT_0106;
    const ccomplex_t IT_0248 = m_N_1*IT_0097;
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = IT_0247 + IT_0249;
    const ccomplex_t IT_0251 = IT_0246*IT_0250;
    const ccomplex_t IT_0252 = IT_0217 + IT_0232 + IT_0234 + -IT_0235 + (-0.5)
      *IT_0236 + IT_0251;
    const ccomplex_t IT_0253 = IT_0038*(IT_0085 + -IT_0086);
    const ccomplex_t IT_0254 = (-0.5)*IT_0253;
    const ccomplex_t IT_0255 = IT_0252 + IT_0254;
    const ccomplex_t IT_0256 = conj(IT_0252) + conj(IT_0254);
    const ccomplex_t IT_0257 = (-24)*IT_0070;
    const ccomplex_t IT_0258 = 6*IT_0068;
    const ccomplex_t IT_0259 = IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = IT_0140*IT_0157;
    const ccomplex_t IT_0261 = IT_0176*IT_0193;
    const ccomplex_t IT_0262 = IT_0035*IT_0055;
    const ccomplex_t IT_0263 = IT_0106*IT_0124;
    const ccomplex_t IT_0264 = IT_0260 + IT_0261 + 0.5*IT_0262 + IT_0263;
    const ccomplex_t IT_0265 = IT_0184*IT_0227;
    const ccomplex_t IT_0266 = IT_0097*IT_0246;
    const ccomplex_t IT_0267 = -IT_0265 + -IT_0266;
    const ccomplex_t IT_0268 = IT_0148*IT_0216;
    const ccomplex_t IT_0269 = IT_0038*IT_0041;
    const ccomplex_t IT_0270 = -IT_0268 + (-0.5)*IT_0269;
    const ccomplex_t IT_0271 = IT_0267 + IT_0270;
    const ccomplex_t IT_0272 = conj(IT_0267) + conj(IT_0270);
    const ccomplex_t IT_0273 = (-12)*IT_0074;
    const ccomplex_t IT_0274 = IT_0071 + IT_0258 + IT_0273;
    const ccomplex_t IT_0275 = 12*IT_0074;
    const ccomplex_t IT_0276 = IT_0069 + IT_0257 + IT_0275;
    const ccomplex_t IT_0277 = m_b*m_N_1;
    const ccomplex_t IT_0278 = pow(s_14, 2);
    const ccomplex_t IT_0279 = IT_0065*IT_0278;
    const ccomplex_t IT_0280 = -IT_0001;
    const ccomplex_t IT_0281 = IT_0279 + IT_0280;
    const ccomplex_t IT_0282 = IT_0277*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0005
      *IT_0006 + 3*IT_0004*IT_0014)*U_sb_00*conj(U_sb_00) + IT_0005*IT_0006
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0284 = (-0.666666666666667)*IT_0283;
    const ccomplex_t IT_0285 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0286 = IT_0053*IT_0285;
    const ccomplex_t IT_0287 = (0 + _Complex_I*1)*IT_0286;
    const ccomplex_t IT_0288 = IT_0284*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0005
      *IT_0006 + 3*IT_0004*IT_0014)*U_sb_00*conj(U_sb_01) + IT_0005*IT_0006
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0290 = (-0.666666666666667)*IT_0289;
    const ccomplex_t IT_0291 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0292 = IT_0005*IT_0291;
    const ccomplex_t IT_0293 = 1.4142135623731*IT_0292;
    const ccomplex_t IT_0294 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0295 = IT_0014*IT_0294;
    const ccomplex_t IT_0296 = 1.4142135623731*IT_0295;
    const ccomplex_t IT_0297 = m_b*conj(N_d1)*e_em*IT_0016*U_sb_11;
    const ccomplex_t IT_0298 = IT_0015*IT_0297;
    const ccomplex_t IT_0299 = 1.4142135623731*IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*(IT_0293 + (-3)*IT_0296 + 3
      *IT_0299);
    const ccomplex_t IT_0301 = 0.166666666666667*IT_0300;
    const ccomplex_t IT_0302 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0303 = IT_0301*IT_0302;
    const ccomplex_t IT_0304 = (0 + _Complex_I*1)*IT_0303;
    const ccomplex_t IT_0305 = IT_0290*IT_0304;
    const ccomplex_t IT_0306 = (-2)*IT_0260 + (-2)*IT_0261 + -IT_0262 + (-2)
      *IT_0263 + -IT_0288 + -IT_0305;
    const ccomplex_t IT_0307 = IT_0021*IT_0285;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*IT_0307;
    const ccomplex_t IT_0309 = IT_0284*IT_0308;
    const ccomplex_t IT_0310 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0311 = IT_0005*IT_0310;
    const ccomplex_t IT_0312 = 1.4142135623731*IT_0311;
    const ccomplex_t IT_0313 = m_b*N_d1*e_em*IT_0016*U_sb_01;
    const ccomplex_t IT_0314 = IT_0015*IT_0313;
    const ccomplex_t IT_0315 = 1.4142135623731*IT_0314;
    const ccomplex_t IT_0316 = (0 + _Complex_I*1)*(IT_0312 + 1.5*IT_0315);
    const ccomplex_t IT_0317 = (-0.333333333333333)*IT_0316;
    const ccomplex_t IT_0318 = IT_0302*IT_0317;
    const ccomplex_t IT_0319 = (0 + _Complex_I*1)*IT_0318;
    const ccomplex_t IT_0320 = IT_0290*IT_0319;
    const ccomplex_t IT_0321 = 2*IT_0265 + 2*IT_0266 + 2*IT_0268 + IT_0269 +
       IT_0309 + IT_0320;
    const ccomplex_t IT_0322 = 6*IT_0321;
    const ccomplex_t IT_0323 = (-6)*IT_0039;
    const ccomplex_t IT_0324 = IT_0322 + IT_0323;
    const ccomplex_t IT_0325 = 6*conj(IT_0321);
    const ccomplex_t IT_0326 = (-6)*conj(IT_0039);
    const ccomplex_t IT_0327 = IT_0325 + IT_0326;
    const ccomplex_t IT_0328 = -IT_0056;
    const ccomplex_t IT_0329 = 6*IT_0039;
    const ccomplex_t IT_0330 = (-6)*IT_0321;
    const ccomplex_t IT_0331 = IT_0329 + IT_0330;
    const ccomplex_t IT_0332 = 6*conj(IT_0039);
    const ccomplex_t IT_0333 = (-6)*conj(IT_0321);
    const ccomplex_t IT_0334 = IT_0332 + IT_0333;
    const ccomplex_t IT_0335 = m_b*IT_0281;
    const ccomplex_t IT_0336 = 6*IT_0306;
    const ccomplex_t IT_0337 = (-6)*IT_0328;
    const ccomplex_t IT_0338 = IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = 6*conj(IT_0306);
    const ccomplex_t IT_0340 = (-6)*conj(IT_0328);
    const ccomplex_t IT_0341 = IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = s_14*s_34*IT_0065;
    const ccomplex_t IT_0343 = -IT_0342;
    const ccomplex_t IT_0344 = s_13 + IT_0343;
    const ccomplex_t IT_0345 = s_13*IT_0344;
    const ccomplex_t IT_0346 = (-2)*IT_0024 + -IT_0309 + -IT_0320;
    const ccomplex_t IT_0347 = 2*IT_0064 + IT_0288 + IT_0305;
    const ccomplex_t IT_0348 = 6*IT_0328;
    const ccomplex_t IT_0349 = (-6)*IT_0306;
    const ccomplex_t IT_0350 = IT_0348 + IT_0349;
    const ccomplex_t IT_0351 = 6*conj(IT_0328);
    const ccomplex_t IT_0352 = (-6)*conj(IT_0306);
    const ccomplex_t IT_0353 = IT_0351 + IT_0352;
    const ccomplex_t IT_0354 = s_13*s_14;
    const ccomplex_t IT_0355 = s_34*IT_0001;
    const ccomplex_t IT_0356 = -IT_0355;
    const ccomplex_t IT_0357 = IT_0354 + IT_0356;
    const ccomplex_t IT_0358 = conj(IT_0025) + conj(IT_0040);
    const ccomplex_t IT_0359 = IT_0025 + IT_0040;
    const ccomplex_t IT_0360 = conj(IT_0057) + conj(IT_0064);
    const ccomplex_t IT_0361 = IT_0057 + IT_0064;
    const ccomplex_t IT_0362 = s_34*IT_0073;
    const ccomplex_t IT_0363 = IT_0065*IT_0362;
    const ccomplex_t IT_0364 = 0.5*IT_0363;
    const ccomplex_t IT_0365 = s_34 + IT_0364;
    const ccomplex_t IT_0366 = m_N_1*IT_0365;
    const ccomplex_t IT_0367 = 2*IT_0366;
    const ccomplex_t IT_0368 = 6*IT_0198;
    const ccomplex_t IT_0369 = 6*IT_0200;
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = 6*conj(IT_0198);
    const ccomplex_t IT_0372 = 6*conj(IT_0200);
    const ccomplex_t IT_0373 = IT_0371 + IT_0372;
    const ccomplex_t IT_0374 = 6*IT_0252;
    const ccomplex_t IT_0375 = 6*IT_0254;
    const ccomplex_t IT_0376 = IT_0374 + IT_0375;
    const ccomplex_t IT_0377 = 6*conj(IT_0252);
    const ccomplex_t IT_0378 = 6*conj(IT_0254);
    const ccomplex_t IT_0379 = IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = s_14*IT_0073;
    const ccomplex_t IT_0381 = IT_0065*IT_0380;
    const ccomplex_t IT_0382 = (-0.25)*IT_0381;
    const ccomplex_t IT_0383 = s_14 + IT_0382;
    const ccomplex_t IT_0384 = m_b*IT_0383;
    const ccomplex_t IT_0385 = (-4)*IT_0384;
    const ccomplex_t IT_0386 = (-0.25)*IT_0363;
    const ccomplex_t IT_0387 = s_34 + IT_0386;
    const ccomplex_t IT_0388 = m_N_1*IT_0387;
    const ccomplex_t IT_0389 = (-4)*IT_0388;
    const ccomplex_t IT_0390 = (-6)*IT_0198;
    const ccomplex_t IT_0391 = (-6)*IT_0200;
    const ccomplex_t IT_0392 = IT_0390 + IT_0391;
    const ccomplex_t IT_0393 = (-6)*conj(IT_0198);
    const ccomplex_t IT_0394 = (-6)*conj(IT_0200);
    const ccomplex_t IT_0395 = IT_0393 + IT_0394;
    const ccomplex_t IT_0396 = (-6)*IT_0252;
    const ccomplex_t IT_0397 = (-6)*IT_0254;
    const ccomplex_t IT_0398 = IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = (-6)*conj(IT_0252);
    const ccomplex_t IT_0400 = (-6)*conj(IT_0254);
    const ccomplex_t IT_0401 = IT_0399 + IT_0400;
    const ccomplex_t IT_0402 = 0.5*IT_0381;
    const ccomplex_t IT_0403 = s_14 + IT_0402;
    const ccomplex_t IT_0404 = m_b*IT_0403;
    const ccomplex_t IT_0405 = 2*IT_0404;
    const ccomplex_t IT_0406 = IT_0073*IT_0277;
    const ccomplex_t IT_0407 = 3*IT_0406;
    const ccomplex_t IT_0408 = 6*IT_0270;
    const ccomplex_t IT_0409 = (-0.166666666666667)*IT_0408;
    const ccomplex_t IT_0410 = 6*IT_0267;
    const ccomplex_t IT_0411 = (-0.166666666666667)*IT_0410;
    const ccomplex_t IT_0412 = IT_0040 + IT_0409 + IT_0411;
    const ccomplex_t IT_0413 = 6*conj(IT_0270);
    const ccomplex_t IT_0414 = (-0.166666666666667)*IT_0413;
    const ccomplex_t IT_0415 = 6*conj(IT_0267);
    const ccomplex_t IT_0416 = (-0.166666666666667)*IT_0415;
    const ccomplex_t IT_0417 = conj(IT_0040) + IT_0414 + IT_0416;
    const ccomplex_t IT_0418 = 6*IT_0025;
    const ccomplex_t IT_0419 = 6*conj(IT_0025);
    const ccomplex_t IT_0420 = (-3)*IT_0406;
    const ccomplex_t IT_0421 = 6*(conj(IT_0057)*(IT_0040 + -IT_0271) + IT_0057
      *(conj(IT_0040) + -IT_0272) + -conj(IT_0264)*IT_0412 + -IT_0264*IT_0417 +
       0.166666666666667*conj(IT_0064)*IT_0418 + 0.166666666666667*IT_0064
      *IT_0419)*IT_0420;
    const ccomplex_t IT_0422 = (-18)*conj(IT_0198);
    const ccomplex_t IT_0423 = s_13*IT_0281;
    const ccomplex_t IT_0424 = IT_0277*IT_0344;
    const ccomplex_t IT_0425 = 6*IT_0346;
    const ccomplex_t IT_0426 = 6*conj(IT_0346);
    const ccomplex_t IT_0427 = (-6)*IT_0346;
    const ccomplex_t IT_0428 = (-6)*conj(IT_0346);
    const ccomplex_t IT_0429 = pow(s_34, 2);
    const ccomplex_t IT_0430 = IT_0065*IT_0429;
    const ccomplex_t IT_0431 = -IT_0430;
    const ccomplex_t IT_0432 = IT_0000 + IT_0431;
    const ccomplex_t IT_0433 = IT_0277*IT_0432;
    const ccomplex_t IT_0434 = s_13*IT_0432;
    const ccomplex_t IT_0435 = (-6)*IT_0347;
    const ccomplex_t IT_0436 = m_b*IT_0344;
    const ccomplex_t IT_0437 = (-6)*conj(IT_0347);
    const ccomplex_t IT_0438 = s_14*IT_0000;
    const ccomplex_t IT_0439 = s_13*s_34;
    const ccomplex_t IT_0440 = -IT_0439;
    const ccomplex_t IT_0441 = IT_0438 + IT_0440;
    const ccomplex_t IT_0442 = conj(IT_0264)*IT_0347 + IT_0264*conj(IT_0347) +
       0.166666666666667*IT_0272*IT_0425 + 0.166666666666667*IT_0271*IT_0426 +
       0.166666666666667*IT_0358*IT_0427 + 0.166666666666667*IT_0359*IT_0428 +
       0.166666666666667*IT_0360*IT_0435 + 0.166666666666667*IT_0361*IT_0437;
    const ccomplex_t IT_0443 = m_N_1*IT_0344;
    const ccomplex_t IT_0444 = m_N_1*IT_0432;
    const ccomplex_t IT_0445 = (conj(IT_0025)*IT_0040 + IT_0025*conj(IT_0040) 
      + conj(IT_0057)*IT_0064 + IT_0057*conj(IT_0064))*IT_0072 + IT_0080*
      (IT_0201*IT_0202 + IT_0255*IT_0256) + IT_0259*(conj(IT_0064)*IT_0264 +
       IT_0064*conj(IT_0264) + conj(IT_0025)*IT_0271 + IT_0025*IT_0272) + 
      (IT_0025*conj(IT_0025) + IT_0040*conj(IT_0040) + IT_0057*conj(IT_0057) +
       IT_0064*conj(IT_0064) + IT_0264*conj(IT_0264) + (IT_0267 + IT_0270)
      *IT_0272)*IT_0274 + (conj(IT_0057)*IT_0264 + IT_0057*conj(IT_0264) + conj
      (IT_0040)*IT_0271 + IT_0040*IT_0272)*IT_0276 + IT_0282*(conj(IT_0306)
      *IT_0324 + IT_0306*IT_0327 + conj(IT_0328)*IT_0331 + IT_0328*IT_0334) +
       IT_0335*(IT_0202*IT_0324 + IT_0201*IT_0327 + IT_0256*IT_0338 + IT_0255
      *IT_0341) + IT_0345*(IT_0334*IT_0346 + IT_0331*conj(IT_0346) + conj
      (IT_0347)*IT_0350 + IT_0347*IT_0353) + IT_0357*(IT_0272*IT_0324 + IT_0271
      *IT_0327 + conj(IT_0264)*IT_0338 + IT_0264*IT_0341 + IT_0331*IT_0358 +
       IT_0334*IT_0359 + IT_0350*IT_0360 + IT_0353*IT_0361) + IT_0367*(conj
      (IT_0064)*IT_0370 + IT_0064*IT_0373 + conj(IT_0025)*IT_0376 + IT_0025
      *IT_0379) + (conj(IT_0025)*IT_0370 + IT_0025*IT_0373 + conj(IT_0064)
      *IT_0376 + IT_0064*IT_0379)*IT_0385 + IT_0389*(conj(IT_0264)*IT_0370 +
       IT_0264*IT_0373 + IT_0272*IT_0376 + IT_0271*IT_0379 + conj(IT_0057)
      *IT_0392 + IT_0057*IT_0395 + conj(IT_0040)*IT_0398 + IT_0040*IT_0401) + 
      (IT_0272*IT_0370 + IT_0271*IT_0373 + conj(IT_0264)*IT_0376 + IT_0264
      *IT_0379 + conj(IT_0040)*IT_0392 + IT_0040*IT_0395 + conj(IT_0057)*IT_0398
       + IT_0057*IT_0401)*IT_0405 + (-6)*IT_0407*(conj(IT_0025)*IT_0057 +
       IT_0025*conj(IT_0057) + conj(IT_0064)*IT_0412 + IT_0064*IT_0417 + (
      -0.166666666666667)*conj(IT_0264)*IT_0418 + (-0.166666666666667)*IT_0264
      *IT_0419) + IT_0421 + IT_0277*(IT_0201*((-18)*conj(IT_0252) + (-18)*conj
      (IT_0254)) + IT_0255*((-18)*conj(IT_0200) + IT_0422)) + (conj(IT_0321)
      *IT_0322 + conj(IT_0039)*IT_0331 + IT_0039*IT_0333 + conj(IT_0306)*IT_0336
       + conj(IT_0328)*IT_0350 + IT_0328*IT_0352)*IT_0423 + IT_0424*(IT_0334
      *IT_0347 + IT_0331*conj(IT_0347) + conj(IT_0328)*IT_0425 + IT_0328*IT_0426
       + conj(IT_0306)*IT_0427 + IT_0306*IT_0428) + (conj(IT_0347)*IT_0427 +
       IT_0347*IT_0428)*IT_0433 + IT_0434*(conj(IT_0346)*IT_0427 + conj(IT_0347)
      *IT_0435) + IT_0436*(IT_0202*IT_0427 + IT_0201*IT_0428 + IT_0256*IT_0435 +
       IT_0255*IT_0437) + 6*IT_0441*IT_0442 + (IT_0256*IT_0331 + IT_0255*IT_0334
       + conj(IT_0200)*IT_0350 + IT_0200*IT_0353 + conj(IT_0328)*IT_0368 +
       IT_0328*IT_0371 + conj(IT_0306)*IT_0390 + IT_0306*IT_0393)*IT_0443 + 
      (conj(IT_0347)*IT_0390 + IT_0347*IT_0393 + IT_0256*IT_0427 + IT_0255
      *IT_0428 + conj(IT_0200)*IT_0435 + IT_0200*IT_0437)*IT_0444;
    return create_ccomplex_return(IT_0445);
}

