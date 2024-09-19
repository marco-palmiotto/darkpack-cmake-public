#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_sb_2_to_b_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_sb_2_to_b_Z(
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = cos(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0006*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = m_b*N_d4*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0016 + 1.5*IT_0023);
    const ccomplex_t IT_0025 = 0.666666666666667*IT_0024;
    const ccomplex_t IT_0026 = 1.5*IT_0025;
    const ccomplex_t IT_0027 = (-0.333333333333333)*IT_0026;
    const ccomplex_t IT_0028 = pow(m_b, 2);
    const ccomplex_t IT_0029 = pow(m_sb_2, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0028 + IT_0029 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0027*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0013*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0012;
    const ccomplex_t IT_0036 = IT_0032*IT_0035;
    const ccomplex_t IT_0037 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = IT_0007*IT_0037;
    const ccomplex_t IT_0040 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0041 = IT_0003*IT_0040;
    const ccomplex_t IT_0042 = IT_0007*IT_0040;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0038 + IT_0039 + 
      -IT_0041 + -IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0046 = IT_0001*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = m_b*N_d1*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0049 = IT_0019*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0047 + 1.5*IT_0050);
    const ccomplex_t IT_0052 = (-0.333333333333333)*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0028 + IT_0029 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0044*IT_0055;
    const ccomplex_t IT_0057 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0058 = IT_0003*IT_0057;
    const ccomplex_t IT_0059 = IT_0007*IT_0057;
    const ccomplex_t IT_0060 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0061 = IT_0007*IT_0060;
    const ccomplex_t IT_0062 = IT_0003*IT_0060;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0058 + IT_0059 + 
      -IT_0061 + -IT_0062);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0069 = IT_0001*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = m_b*N_d2*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0072 = IT_0019*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0070 + 1.5*IT_0073);
    const ccomplex_t IT_0075 = 0.666666666666667*IT_0074;
    const ccomplex_t IT_0076 = 1.5*IT_0075;
    const ccomplex_t IT_0077 = (-0.333333333333333)*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0028 + IT_0029 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0067*IT_0080;
    const ccomplex_t IT_0082 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0083 = IT_0003*IT_0082;
    const ccomplex_t IT_0084 = IT_0007*IT_0082;
    const ccomplex_t IT_0085 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0086 = IT_0003*IT_0085;
    const ccomplex_t IT_0087 = IT_0007*IT_0085;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0083 + IT_0084 + 
      -IT_0086 + -IT_0087);
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0091 = IT_0001*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = m_b*N_d3*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0094 = IT_0019*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0092 + 1.5*IT_0095);
    const ccomplex_t IT_0097 = 0.666666666666667*IT_0096;
    const ccomplex_t IT_0098 = 1.5*IT_0097;
    const ccomplex_t IT_0099 = 0.666666666666667*IT_0098;
    const ccomplex_t IT_0100 = 1.5*IT_0099;
    const ccomplex_t IT_0101 = 0.666666666666667*IT_0100;
    const ccomplex_t IT_0102 = 1.5*IT_0101;
    const ccomplex_t IT_0103 = (-0.333333333333333)*IT_0102;
    const ccomplex_t IT_0104 = cpow((-2)*s_23 + IT_0028 + IT_0029 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0089*IT_0106;
    const ccomplex_t IT_0108 = (-0.5)*IT_0036 + -IT_0056 + -IT_0081 + -IT_0107;
    const ccomplex_t IT_0109 = m_b*conj(N_d4)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0110 = IT_0019*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0113 = IT_0001*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0116 = IT_0006*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0111 + 0.333333333333333
      *IT_0114 + -IT_0117);
    const ccomplex_t IT_0119 = 3*IT_0118;
    const ccomplex_t IT_0120 = 0.166666666666667*IT_0119;
    const ccomplex_t IT_0121 = IT_0030*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = IT_0035*IT_0122;
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = IT_0013*IT_0122;
    const ccomplex_t IT_0126 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0127 = IT_0003*IT_0126;
    const ccomplex_t IT_0128 = IT_0007*IT_0126;
    const ccomplex_t IT_0129 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0130 = IT_0003*IT_0129;
    const ccomplex_t IT_0131 = IT_0007*IT_0129;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0127 + IT_0128 + 
      -IT_0130 + -IT_0131);
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = m_b*conj(N_d1)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0135 = IT_0019*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0138 = IT_0006*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0141 = IT_0001*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0136 + -IT_0139 +
       0.333333333333333*IT_0142);
    const ccomplex_t IT_0144 = 3*IT_0143;
    const ccomplex_t IT_0145 = 0.166666666666667*IT_0144;
    const ccomplex_t IT_0146 = IT_0053*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = IT_0133*IT_0147;
    const ccomplex_t IT_0149 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0150 = IT_0003*IT_0149;
    const ccomplex_t IT_0151 = IT_0007*IT_0149;
    const ccomplex_t IT_0152 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0153 = IT_0007*IT_0152;
    const ccomplex_t IT_0154 = IT_0003*IT_0152;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0150 + IT_0151 + 
      -IT_0153 + -IT_0154);
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = (-0.5)*IT_0156;
    const ccomplex_t IT_0158 = m_b*conj(N_d2)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0159 = IT_0019*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0162 = IT_0001*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0165 = IT_0006*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0160 + 0.333333333333333
      *IT_0163 + -IT_0166);
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = (-3)*IT_0168;
    const ccomplex_t IT_0170 = 0.166666666666667*IT_0169;
    const ccomplex_t IT_0171 = IT_0078*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = IT_0157*IT_0172;
    const ccomplex_t IT_0174 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0175 = IT_0003*IT_0174;
    const ccomplex_t IT_0176 = IT_0007*IT_0174;
    const ccomplex_t IT_0177 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0178 = IT_0003*IT_0177;
    const ccomplex_t IT_0179 = IT_0007*IT_0177;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*(IT_0175 + IT_0176 + 
      -IT_0178 + -IT_0179);
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = m_b*conj(N_d3)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0183 = IT_0019*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0186 = IT_0001*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0189 = IT_0006*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*(IT_0184 + 0.333333333333333
      *IT_0187 + -IT_0190);
    const ccomplex_t IT_0192 = 3*IT_0191;
    const ccomplex_t IT_0193 = 0.166666666666667*IT_0192;
    const ccomplex_t IT_0194 = IT_0104*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = IT_0181*IT_0195;
    const ccomplex_t IT_0197 = 0.5*IT_0125 + IT_0148 + IT_0173 + IT_0196;
    const ccomplex_t IT_0198 = pow(m_Z, -2);
    const ccomplex_t IT_0199 = pow(m_Z, 4);
    const ccomplex_t IT_0200 = s_13*IT_0199;
    const ccomplex_t IT_0201 = IT_0198*IT_0200;
    const ccomplex_t IT_0202 = (-6)*IT_0201;
    const ccomplex_t IT_0203 = s_14*s_34;
    const ccomplex_t IT_0204 = (-24)*IT_0203;
    const ccomplex_t IT_0205 = pow(m_Z, 2);
    const ccomplex_t IT_0206 = s_13*IT_0205;
    const ccomplex_t IT_0207 = 12*IT_0206;
    const ccomplex_t IT_0208 = IT_0202 + IT_0204 + IT_0207;
    const ccomplex_t IT_0209 = 6*IT_0201;
    const ccomplex_t IT_0210 = 24*IT_0203;
    const ccomplex_t IT_0211 = (-12)*IT_0206;
    const ccomplex_t IT_0212 = IT_0209 + IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = (-2)*IT_0203;
    const ccomplex_t IT_0214 = IT_0206 + IT_0213;
    const ccomplex_t IT_0215 = IT_0198*IT_0214;
    const ccomplex_t IT_0216 = (-6)*IT_0215;
    const ccomplex_t IT_0217 = 12*s_13;
    const ccomplex_t IT_0218 = IT_0216 + IT_0217;
    const ccomplex_t IT_0219 = pow(m_N_4, 2);
    const ccomplex_t IT_0220 = cpow((-2)*s_12 + IT_0028 + -IT_0029 + -IT_0219 
      + -reg_prop, -1);
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0222 = e_em*IT_0003;
    const ccomplex_t IT_0223 = e_em*IT_0007;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*(IT_0222 + 3*IT_0223);
    const ccomplex_t IT_0225 = (-0.166666666666667)*IT_0224;
    const ccomplex_t IT_0226 = IT_0027*IT_0225;
    const ccomplex_t IT_0227 = IT_0220*IT_0221*IT_0226;
    const ccomplex_t IT_0228 = m_N_4*IT_0013;
    const ccomplex_t IT_0229 = m_N_4*IT_0035;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = IT_0122*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0234 = 0.333333333333333*IT_0233;
    const ccomplex_t IT_0235 = IT_0027*IT_0234;
    const ccomplex_t IT_0236 = IT_0220*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = m_b*IT_0237;
    const ccomplex_t IT_0239 = -IT_0227 + (-0.5)*IT_0232 + IT_0238;
    const ccomplex_t IT_0240 = m_N_1*IT_0044;
    const ccomplex_t IT_0241 = m_N_4*IT_0133;
    const ccomplex_t IT_0242 = -IT_0241;
    const ccomplex_t IT_0243 = IT_0240 + IT_0242;
    const ccomplex_t IT_0244 = IT_0147*IT_0243;
    const ccomplex_t IT_0245 = m_N_2*IT_0067;
    const ccomplex_t IT_0246 = m_N_4*IT_0157;
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = IT_0245 + IT_0247;
    const ccomplex_t IT_0249 = IT_0172*IT_0248;
    const ccomplex_t IT_0250 = m_N_3*IT_0089;
    const ccomplex_t IT_0251 = m_N_4*IT_0181;
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = IT_0250 + IT_0252;
    const ccomplex_t IT_0254 = IT_0195*IT_0253;
    const ccomplex_t IT_0255 = (-0.5)*IT_0232 + IT_0244 + IT_0249 + IT_0254;
    const ccomplex_t IT_0256 = IT_0239 + IT_0255;
    const ccomplex_t IT_0257 = conj(IT_0239) + conj(IT_0255);
    const ccomplex_t IT_0258 = IT_0120*IT_0234;
    const ccomplex_t IT_0259 = IT_0220*IT_0221*IT_0258;
    const ccomplex_t IT_0260 = IT_0120*IT_0225;
    const ccomplex_t IT_0261 = IT_0220*IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = m_b*IT_0262;
    const ccomplex_t IT_0264 = IT_0032*IT_0231;
    const ccomplex_t IT_0265 = IT_0259 + -IT_0263 + (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = m_N_1*IT_0133;
    const ccomplex_t IT_0267 = m_N_4*IT_0044;
    const ccomplex_t IT_0268 = -IT_0267;
    const ccomplex_t IT_0269 = IT_0266 + IT_0268;
    const ccomplex_t IT_0270 = IT_0055*IT_0269;
    const ccomplex_t IT_0271 = m_N_2*IT_0157;
    const ccomplex_t IT_0272 = m_N_4*IT_0067;
    const ccomplex_t IT_0273 = -IT_0272;
    const ccomplex_t IT_0274 = IT_0271 + IT_0273;
    const ccomplex_t IT_0275 = IT_0080*IT_0274;
    const ccomplex_t IT_0276 = m_N_3*IT_0181;
    const ccomplex_t IT_0277 = m_N_4*IT_0089;
    const ccomplex_t IT_0278 = -IT_0277;
    const ccomplex_t IT_0279 = IT_0276 + IT_0278;
    const ccomplex_t IT_0280 = IT_0106*IT_0279;
    const ccomplex_t IT_0281 = (-0.5)*IT_0264 + -IT_0270 + -IT_0275 + -IT_0280;
    const ccomplex_t IT_0282 = IT_0265 + IT_0281;
    const ccomplex_t IT_0283 = conj(IT_0265) + conj(IT_0281);
    const ccomplex_t IT_0284 = m_b*m_N_4;
    const ccomplex_t IT_0285 = pow(s_14, 2);
    const ccomplex_t IT_0286 = IT_0198*IT_0285;
    const ccomplex_t IT_0287 = -IT_0219;
    const ccomplex_t IT_0288 = IT_0286 + IT_0287;
    const ccomplex_t IT_0289 = IT_0284*IT_0288;
    const ccomplex_t IT_0290 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0006)*conj(U_sb_00)*U_sb_01 + IT_0001*IT_0002*conj
      (U_sb_10)*U_sb_11);
    const ccomplex_t IT_0291 = (-0.666666666666667)*IT_0290;
    const ccomplex_t IT_0292 = m_b*conj(N_d4)*e_em*IT_0020*U_sb_10;
    const ccomplex_t IT_0293 = IT_0019*IT_0292;
    const ccomplex_t IT_0294 = 1.4142135623731*IT_0293;
    const ccomplex_t IT_0295 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0296 = IT_0006*IT_0295;
    const ccomplex_t IT_0297 = 1.4142135623731*IT_0296;
    const ccomplex_t IT_0298 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0299 = IT_0001*IT_0298;
    const ccomplex_t IT_0300 = 1.4142135623731*IT_0299;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*(IT_0294 + -IT_0297 +
       0.333333333333333*IT_0300);
    const ccomplex_t IT_0302 = 3*IT_0301;
    const ccomplex_t IT_0303 = 0.333333333333333*IT_0302;
    const ccomplex_t IT_0304 = 3*IT_0303;
    const ccomplex_t IT_0305 = 0.166666666666667*IT_0304;
    const ccomplex_t IT_0306 = cpow((-2)*s_13 + IT_0028 + IT_0219 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0307 = IT_0305*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*IT_0307;
    const ccomplex_t IT_0309 = IT_0291*IT_0308;
    const ccomplex_t IT_0310 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0006)*U_sb_01*conj(U_sb_01) + IT_0001*IT_0002
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0311 = (-0.666666666666667)*IT_0310;
    const ccomplex_t IT_0312 = cpow((-2)*s_13 + IT_0028 + IT_0219 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0313 = IT_0120*IT_0312;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*IT_0313;
    const ccomplex_t IT_0315 = IT_0311*IT_0314;
    const ccomplex_t IT_0316 = -IT_0123 + -IT_0309 + -IT_0315;
    const ccomplex_t IT_0317 = IT_0027*IT_0312;
    const ccomplex_t IT_0318 = (0 + _Complex_I*1)*IT_0317;
    const ccomplex_t IT_0319 = IT_0311*IT_0318;
    const ccomplex_t IT_0320 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0321 = IT_0001*IT_0320;
    const ccomplex_t IT_0322 = 1.4142135623731*IT_0321;
    const ccomplex_t IT_0323 = m_b*N_d4*e_em*IT_0020*U_sb_00;
    const ccomplex_t IT_0324 = IT_0019*IT_0323;
    const ccomplex_t IT_0325 = 1.4142135623731*IT_0324;
    const ccomplex_t IT_0326 = (0 + _Complex_I*1)*(IT_0322 + 1.5*IT_0325);
    const ccomplex_t IT_0327 = 0.666666666666667*IT_0326;
    const ccomplex_t IT_0328 = 1.5*IT_0327;
    const ccomplex_t IT_0329 = (-0.333333333333333)*IT_0328;
    const ccomplex_t IT_0330 = IT_0306*IT_0329;
    const ccomplex_t IT_0331 = (0 + _Complex_I*1)*IT_0330;
    const ccomplex_t IT_0332 = IT_0291*IT_0331;
    const ccomplex_t IT_0333 = IT_0033 + IT_0319 + IT_0332;
    const ccomplex_t IT_0334 = 6*IT_0333;
    const ccomplex_t IT_0335 = IT_0036 + 2*IT_0056 + 2*IT_0081 + 2*IT_0107;
    const ccomplex_t IT_0336 = (-6)*IT_0335;
    const ccomplex_t IT_0337 = IT_0334 + IT_0336;
    const ccomplex_t IT_0338 = 6*conj(IT_0333);
    const ccomplex_t IT_0339 = (-6)*conj(IT_0335);
    const ccomplex_t IT_0340 = IT_0338 + IT_0339;
    const ccomplex_t IT_0341 = -IT_0125 + (-2)*IT_0148 + (-2)*IT_0173 + (-2)
      *IT_0196;
    const ccomplex_t IT_0342 = 6*IT_0335;
    const ccomplex_t IT_0343 = (-6)*IT_0333;
    const ccomplex_t IT_0344 = IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = 6*conj(IT_0335);
    const ccomplex_t IT_0346 = (-6)*conj(IT_0333);
    const ccomplex_t IT_0347 = IT_0345 + IT_0346;
    const ccomplex_t IT_0348 = m_b*IT_0288;
    const ccomplex_t IT_0349 = 6*IT_0316;
    const ccomplex_t IT_0350 = (-6)*IT_0341;
    const ccomplex_t IT_0351 = IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = 6*conj(IT_0316);
    const ccomplex_t IT_0353 = (-6)*conj(IT_0341);
    const ccomplex_t IT_0354 = IT_0352 + IT_0353;
    const ccomplex_t IT_0355 = s_14*s_34*IT_0198;
    const ccomplex_t IT_0356 = -IT_0355;
    const ccomplex_t IT_0357 = s_13 + IT_0356;
    const ccomplex_t IT_0358 = s_13*IT_0357;
    const ccomplex_t IT_0359 = (-2)*IT_0237 + -IT_0319 + -IT_0332;
    const ccomplex_t IT_0360 = 2*IT_0262 + IT_0309 + IT_0315;
    const ccomplex_t IT_0361 = 6*IT_0341;
    const ccomplex_t IT_0362 = (-6)*IT_0316;
    const ccomplex_t IT_0363 = IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = 6*conj(IT_0341);
    const ccomplex_t IT_0365 = (-6)*conj(IT_0316);
    const ccomplex_t IT_0366 = IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = s_13*s_14;
    const ccomplex_t IT_0368 = s_34*IT_0219;
    const ccomplex_t IT_0369 = -IT_0368;
    const ccomplex_t IT_0370 = IT_0367 + IT_0369;
    const ccomplex_t IT_0371 = -IT_0237;
    const ccomplex_t IT_0372 = conj(IT_0108) + conj(IT_0371);
    const ccomplex_t IT_0373 = IT_0108 + IT_0371;
    const ccomplex_t IT_0374 = conj(IT_0197) + conj(IT_0262);
    const ccomplex_t IT_0375 = IT_0197 + IT_0262;
    const ccomplex_t IT_0376 = s_34*IT_0205;
    const ccomplex_t IT_0377 = IT_0198*IT_0376;
    const ccomplex_t IT_0378 = (-0.25)*IT_0377;
    const ccomplex_t IT_0379 = s_34 + IT_0378;
    const ccomplex_t IT_0380 = m_N_4*IT_0379;
    const ccomplex_t IT_0381 = (-4)*IT_0380;
    const ccomplex_t IT_0382 = 6*IT_0239;
    const ccomplex_t IT_0383 = 6*IT_0255;
    const ccomplex_t IT_0384 = IT_0382 + IT_0383;
    const ccomplex_t IT_0385 = 6*conj(IT_0239);
    const ccomplex_t IT_0386 = 6*conj(IT_0255);
    const ccomplex_t IT_0387 = IT_0385 + IT_0386;
    const ccomplex_t IT_0388 = 6*IT_0265;
    const ccomplex_t IT_0389 = 6*IT_0281;
    const ccomplex_t IT_0390 = IT_0388 + IT_0389;
    const ccomplex_t IT_0391 = 6*conj(IT_0265);
    const ccomplex_t IT_0392 = 6*conj(IT_0281);
    const ccomplex_t IT_0393 = IT_0391 + IT_0392;
    const ccomplex_t IT_0394 = (-6)*IT_0239;
    const ccomplex_t IT_0395 = (-6)*IT_0255;
    const ccomplex_t IT_0396 = IT_0394 + IT_0395;
    const ccomplex_t IT_0397 = (-6)*conj(IT_0239);
    const ccomplex_t IT_0398 = (-6)*conj(IT_0255);
    const ccomplex_t IT_0399 = IT_0397 + IT_0398;
    const ccomplex_t IT_0400 = (-6)*IT_0265;
    const ccomplex_t IT_0401 = (-6)*IT_0281;
    const ccomplex_t IT_0402 = IT_0400 + IT_0401;
    const ccomplex_t IT_0403 = (-6)*conj(IT_0265);
    const ccomplex_t IT_0404 = (-6)*conj(IT_0281);
    const ccomplex_t IT_0405 = IT_0403 + IT_0404;
    const ccomplex_t IT_0406 = s_14*IT_0205;
    const ccomplex_t IT_0407 = IT_0198*IT_0406;
    const ccomplex_t IT_0408 = 0.5*IT_0407;
    const ccomplex_t IT_0409 = s_14 + IT_0408;
    const ccomplex_t IT_0410 = m_b*IT_0409;
    const ccomplex_t IT_0411 = 2*IT_0410;
    const ccomplex_t IT_0412 = IT_0205*IT_0284;
    const ccomplex_t IT_0413 = 3*IT_0412;
    const ccomplex_t IT_0414 = 6*IT_0371;
    const ccomplex_t IT_0415 = 6*conj(IT_0371);
    const ccomplex_t IT_0416 = (-18)*conj(IT_0239);
    const ccomplex_t IT_0417 = s_13*IT_0288;
    const ccomplex_t IT_0418 = IT_0284*IT_0357;
    const ccomplex_t IT_0419 = 6*IT_0359;
    const ccomplex_t IT_0420 = 6*conj(IT_0359);
    const ccomplex_t IT_0421 = (-6)*IT_0359;
    const ccomplex_t IT_0422 = (-6)*conj(IT_0359);
    const ccomplex_t IT_0423 = pow(s_34, 2);
    const ccomplex_t IT_0424 = IT_0198*IT_0423;
    const ccomplex_t IT_0425 = -IT_0424;
    const ccomplex_t IT_0426 = IT_0028 + IT_0425;
    const ccomplex_t IT_0427 = IT_0284*IT_0426;
    const ccomplex_t IT_0428 = s_13*IT_0426;
    const ccomplex_t IT_0429 = (-6)*IT_0360;
    const ccomplex_t IT_0430 = m_b*IT_0357;
    const ccomplex_t IT_0431 = (-6)*conj(IT_0360);
    const ccomplex_t IT_0432 = s_14*IT_0028;
    const ccomplex_t IT_0433 = s_13*s_34;
    const ccomplex_t IT_0434 = -IT_0433;
    const ccomplex_t IT_0435 = IT_0432 + IT_0434;
    const ccomplex_t IT_0436 = conj(IT_0124)*IT_0360 + IT_0124*conj(IT_0360) +
       0.166666666666667*conj(IT_0034)*IT_0419 + 0.166666666666667*IT_0034
      *IT_0420 + 0.166666666666667*IT_0372*IT_0421 + 0.166666666666667*IT_0373
      *IT_0422 + 0.166666666666667*IT_0374*IT_0429 + 0.166666666666667*IT_0375
      *IT_0431;
    const ccomplex_t IT_0437 = m_N_4*IT_0357;
    const ccomplex_t IT_0438 = m_N_4*IT_0426;
    const ccomplex_t IT_0439 = (-3)*IT_0412;
    const ccomplex_t IT_0440 = 6*IT_0034;
    const ccomplex_t IT_0441 = (-0.166666666666667)*IT_0440;
    const ccomplex_t IT_0442 = IT_0108 + IT_0441;
    const ccomplex_t IT_0443 = 6*conj(IT_0034);
    const ccomplex_t IT_0444 = (-0.166666666666667)*IT_0443;
    const ccomplex_t IT_0445 = conj(IT_0108) + IT_0444;
    const ccomplex_t IT_0446 = (-0.166666666666667)*conj(IT_0262);
    const ccomplex_t IT_0447 = IT_0202 + IT_0210;
    const ccomplex_t IT_0448 = IT_0204 + IT_0209;
    const ccomplex_t IT_0449 = (-0.25)*IT_0407;
    const ccomplex_t IT_0450 = s_14 + IT_0449;
    const ccomplex_t IT_0451 = m_b*IT_0450;
    const ccomplex_t IT_0452 = (-4)*IT_0451;
    const ccomplex_t IT_0453 = (-0.166666666666667)*IT_0452;
    const ccomplex_t IT_0454 = 0.5*IT_0377;
    const ccomplex_t IT_0455 = s_34 + IT_0454;
    const ccomplex_t IT_0456 = m_N_4*IT_0455;
    const ccomplex_t IT_0457 = 2*IT_0456;
    const ccomplex_t IT_0458 = (-0.166666666666667)*IT_0457;
    const ccomplex_t IT_0459 = (conj(IT_0034)*IT_0108 + IT_0034*conj(IT_0108) 
      + conj(IT_0124)*IT_0197 + IT_0124*conj(IT_0197))*IT_0208 + (IT_0034*conj
      (IT_0034) + IT_0108*conj(IT_0108) + IT_0124*conj(IT_0124) + IT_0197*conj
      (IT_0197))*IT_0212 + IT_0218*(IT_0256*IT_0257 + IT_0282*IT_0283) + IT_0289
      *(conj(IT_0316)*IT_0337 + IT_0316*IT_0340 + conj(IT_0341)*IT_0344 +
       IT_0341*IT_0347) + IT_0348*(IT_0257*IT_0337 + IT_0256*IT_0340 + IT_0283
      *IT_0351 + IT_0282*IT_0354) + IT_0358*(IT_0347*IT_0359 + IT_0344*conj
      (IT_0359) + conj(IT_0360)*IT_0363 + IT_0360*IT_0366) + IT_0370*(conj
      (IT_0034)*IT_0337 + IT_0034*IT_0340 + conj(IT_0124)*IT_0351 + IT_0124
      *IT_0354 + IT_0344*IT_0372 + IT_0347*IT_0373 + IT_0363*IT_0374 + IT_0366
      *IT_0375) + IT_0381*(conj(IT_0124)*IT_0384 + IT_0124*IT_0387 + conj
      (IT_0034)*IT_0390 + IT_0034*IT_0393 + conj(IT_0197)*IT_0396 + IT_0197
      *IT_0399 + conj(IT_0108)*IT_0402 + IT_0108*IT_0405) + (conj(IT_0034)
      *IT_0384 + IT_0034*IT_0387 + conj(IT_0124)*IT_0390 + IT_0124*IT_0393 +
       conj(IT_0108)*IT_0396 + IT_0108*IT_0399 + conj(IT_0197)*IT_0402 + IT_0197
      *IT_0405)*IT_0411 + IT_0413*(conj(IT_0124)*IT_0414 + IT_0124*IT_0415) +
       IT_0284*(IT_0256*((-18)*conj(IT_0265) + (-18)*conj(IT_0281)) + IT_0282*((
      -18)*conj(IT_0255) + IT_0416)) + (conj(IT_0333)*IT_0334 + conj(IT_0335)
      *IT_0344 + IT_0335*IT_0346 + conj(IT_0316)*IT_0349 + conj(IT_0341)*IT_0363
       + IT_0341*IT_0365)*IT_0417 + IT_0418*(IT_0347*IT_0360 + IT_0344*conj
      (IT_0360) + conj(IT_0341)*IT_0419 + IT_0341*IT_0420 + conj(IT_0316)
      *IT_0421 + IT_0316*IT_0422) + (conj(IT_0360)*IT_0421 + IT_0360*IT_0422)
      *IT_0427 + IT_0428*(conj(IT_0359)*IT_0421 + conj(IT_0360)*IT_0429) +
       IT_0430*(IT_0257*IT_0421 + IT_0256*IT_0422 + IT_0283*IT_0429 + IT_0282
      *IT_0431) + 6*IT_0435*IT_0436 + (IT_0283*IT_0344 + IT_0282*IT_0347 + conj
      (IT_0255)*IT_0363 + IT_0255*IT_0366 + conj(IT_0341)*IT_0382 + IT_0341
      *IT_0385 + conj(IT_0316)*IT_0394 + IT_0316*IT_0397)*IT_0437 + (conj
      (IT_0360)*IT_0394 + IT_0360*IT_0397 + IT_0283*IT_0421 + IT_0282*IT_0422 +
       conj(IT_0255)*IT_0429 + IT_0255*IT_0431)*IT_0438 + (-6)*IT_0439*((conj
      (IT_0034) + -conj(IT_0108))*IT_0197 + (IT_0034 + -IT_0108)*conj(IT_0197) +
       (-0.166666666666667)*IT_0262*IT_0415 + conj(IT_0124)*IT_0442 + IT_0124
      *IT_0445 + IT_0414*IT_0446) + (-6)*conj(IT_0262)*(IT_0413*IT_0442 + -1./6
      *IT_0197*IT_0447 + -1./6*IT_0124*IT_0448 + IT_0390*IT_0453 + IT_0384
      *IT_0458) + (-6)*IT_0262*(IT_0413*IT_0445 + IT_0212*IT_0446 + -1./6*conj
      (IT_0197)*IT_0447 + -1./6*conj(IT_0124)*IT_0448 + IT_0393*IT_0453 +
       IT_0387*IT_0458) + (-6)*conj(IT_0371)*((-0.166666666666667)*IT_0212
      *IT_0371 + IT_0197*IT_0413 + (-0.166666666666667)*IT_0108*IT_0447 + (
      -0.166666666666667)*IT_0034*IT_0448 + IT_0384*IT_0453 + IT_0390*IT_0458) +
       (-6)*IT_0371*(conj(IT_0197)*IT_0413 + (-0.166666666666667)*conj(IT_0108)
      *IT_0447 + (-0.166666666666667)*conj(IT_0034)*IT_0448 + IT_0387*IT_0453 +
       IT_0393*IT_0458);
    return create_ccomplex_return(IT_0459);
}

