#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_sb_1_to_b_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_sb_1_to_b_Z(
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
    const ccomplex_t IT_0004 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = cos(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0006*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = m_b*N_d4*e_em*IT_0020*U_sb_00;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0016 + 1.5*IT_0023);
    const ccomplex_t IT_0025 = (-0.333333333333333)*IT_0024;
    const ccomplex_t IT_0026 = pow(m_b, 2);
    const ccomplex_t IT_0027 = pow(m_sb_1, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0013*IT_0030;
    const ccomplex_t IT_0032 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = IT_0007*IT_0032;
    const ccomplex_t IT_0035 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0036 = IT_0003*IT_0035;
    const ccomplex_t IT_0037 = IT_0007*IT_0035;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0033 + IT_0034 + 
      -IT_0036 + -IT_0037);
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = m_b*N_d2*e_em*IT_0020*U_sb_00;
    const ccomplex_t IT_0044 = IT_0019*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 1.5*IT_0045);
    const ccomplex_t IT_0047 = 0.666666666666667*IT_0046;
    const ccomplex_t IT_0048 = 1.5*IT_0047;
    const ccomplex_t IT_0049 = 0.666666666666667*IT_0048;
    const ccomplex_t IT_0050 = 1.5*IT_0049;
    const ccomplex_t IT_0051 = 0.666666666666667*IT_0050;
    const ccomplex_t IT_0052 = 1.5*IT_0051;
    const ccomplex_t IT_0053 = 0.666666666666667*IT_0052;
    const ccomplex_t IT_0054 = 1.5*IT_0053;
    const ccomplex_t IT_0055 = (-0.333333333333333)*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0039*IT_0058;
    const ccomplex_t IT_0060 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0061 = IT_0003*IT_0060;
    const ccomplex_t IT_0062 = IT_0007*IT_0060;
    const ccomplex_t IT_0063 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0064 = IT_0003*IT_0063;
    const ccomplex_t IT_0065 = IT_0007*IT_0063;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0061 + IT_0062 + 
      -IT_0064 + -IT_0065);
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0069 = IT_0001*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = m_b*N_d3*e_em*IT_0020*U_sb_00;
    const ccomplex_t IT_0072 = IT_0019*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0070 + 1.5*IT_0073);
    const ccomplex_t IT_0075 = 0.666666666666667*IT_0074;
    const ccomplex_t IT_0076 = 1.5*IT_0075;
    const ccomplex_t IT_0077 = (-0.333333333333333)*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0067*IT_0080;
    const ccomplex_t IT_0082 = -IT_0031 + (-0.5)*IT_0059 + -IT_0081;
    const ccomplex_t IT_0083 = 0.5*IT_0038;
    const ccomplex_t IT_0084 = IT_0058*IT_0083;
    const ccomplex_t IT_0085 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0086 = IT_0003*IT_0085;
    const ccomplex_t IT_0087 = IT_0007*IT_0085;
    const ccomplex_t IT_0088 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0089 = IT_0003*IT_0088;
    const ccomplex_t IT_0090 = IT_0007*IT_0088;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0086 + IT_0087 + 
      -IT_0089 + -IT_0090);
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0094 = IT_0001*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = m_b*N_d1*e_em*IT_0020*U_sb_00;
    const ccomplex_t IT_0097 = IT_0019*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0095 + 1.5*IT_0098);
    const ccomplex_t IT_0100 = (-0.333333333333333)*IT_0099;
    const ccomplex_t IT_0101 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0092*IT_0103;
    const ccomplex_t IT_0105 = (-0.5)*IT_0084 + -IT_0104;
    const ccomplex_t IT_0106 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0107 = IT_0003*IT_0106;
    const ccomplex_t IT_0108 = IT_0007*IT_0106;
    const ccomplex_t IT_0109 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0110 = IT_0003*IT_0109;
    const ccomplex_t IT_0111 = IT_0007*IT_0109;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0107 + IT_0108 + 
      -IT_0110 + -IT_0111);
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = m_b*conj(N_d4)*e_em*IT_0020*U_sb_10;
    const ccomplex_t IT_0115 = IT_0019*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0118 = IT_0006*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0121 = IT_0001*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0116 + -IT_0119 +
       0.333333333333333*IT_0122);
    const ccomplex_t IT_0124 = 3*IT_0123;
    const ccomplex_t IT_0125 = 0.166666666666667*IT_0124;
    const ccomplex_t IT_0126 = IT_0028*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = IT_0113*IT_0127;
    const ccomplex_t IT_0129 = m_b*conj(N_d2)*e_em*IT_0020*U_sb_10;
    const ccomplex_t IT_0130 = IT_0019*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0133 = IT_0001*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0136 = IT_0006*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0131 + 0.333333333333333
      *IT_0134 + -IT_0137);
    const ccomplex_t IT_0139 = 3*IT_0138;
    const ccomplex_t IT_0140 = 0.333333333333333*IT_0139;
    const ccomplex_t IT_0141 = 3*IT_0140;
    const ccomplex_t IT_0142 = (-0.333333333333333)*IT_0141;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = 3*IT_0143;
    const ccomplex_t IT_0145 = 0.166666666666667*IT_0144;
    const ccomplex_t IT_0146 = IT_0056*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = IT_0083*IT_0147;
    const ccomplex_t IT_0149 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0150 = IT_0003*IT_0149;
    const ccomplex_t IT_0151 = IT_0007*IT_0149;
    const ccomplex_t IT_0152 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0153 = IT_0003*IT_0152;
    const ccomplex_t IT_0154 = IT_0007*IT_0152;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0150 + IT_0151 + 
      -IT_0153 + -IT_0154);
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = m_b*conj(N_d3)*e_em*IT_0020*U_sb_10;
    const ccomplex_t IT_0158 = IT_0019*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0161 = IT_0006*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0164 = IT_0001*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*(IT_0159 + -IT_0162 +
       0.333333333333333*IT_0165);
    const ccomplex_t IT_0167 = 3*IT_0166;
    const ccomplex_t IT_0168 = 0.166666666666667*IT_0167;
    const ccomplex_t IT_0169 = IT_0078*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = IT_0156*IT_0170;
    const ccomplex_t IT_0172 = IT_0128 + 0.5*IT_0148 + IT_0171;
    const ccomplex_t IT_0173 = IT_0039*IT_0147;
    const ccomplex_t IT_0174 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0175 = IT_0003*IT_0174;
    const ccomplex_t IT_0176 = IT_0007*IT_0174;
    const ccomplex_t IT_0177 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0178 = IT_0003*IT_0177;
    const ccomplex_t IT_0179 = IT_0007*IT_0177;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*(IT_0175 + IT_0176 + 
      -IT_0178 + -IT_0179);
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = m_b*conj(N_d1)*e_em*IT_0020*U_sb_10;
    const ccomplex_t IT_0183 = IT_0019*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0186 = IT_0001*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0189 = IT_0006*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*(IT_0184 + 0.333333333333333
      *IT_0187 + -IT_0190);
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = (-3)*IT_0192;
    const ccomplex_t IT_0194 = (-0.333333333333333)*IT_0193;
    const ccomplex_t IT_0195 = -IT_0194;
    const ccomplex_t IT_0196 = 3*IT_0195;
    const ccomplex_t IT_0197 = 0.166666666666667*IT_0196;
    const ccomplex_t IT_0198 = IT_0101*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = IT_0181*IT_0199;
    const ccomplex_t IT_0201 = 0.5*IT_0173 + IT_0200;
    const ccomplex_t IT_0202 = pow(m_Z, -2);
    const ccomplex_t IT_0203 = pow(m_Z, 4);
    const ccomplex_t IT_0204 = s_13*IT_0203;
    const ccomplex_t IT_0205 = IT_0202*IT_0204;
    const ccomplex_t IT_0206 = (-6)*IT_0205;
    const ccomplex_t IT_0207 = s_14*s_34;
    const ccomplex_t IT_0208 = (-24)*IT_0207;
    const ccomplex_t IT_0209 = pow(m_Z, 2);
    const ccomplex_t IT_0210 = s_13*IT_0209;
    const ccomplex_t IT_0211 = 12*IT_0210;
    const ccomplex_t IT_0212 = IT_0206 + IT_0208 + IT_0211;
    const ccomplex_t IT_0213 = 6*IT_0205;
    const ccomplex_t IT_0214 = 24*IT_0207;
    const ccomplex_t IT_0215 = (-12)*IT_0210;
    const ccomplex_t IT_0216 = IT_0213 + IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = (-2)*IT_0207;
    const ccomplex_t IT_0218 = IT_0210 + IT_0217;
    const ccomplex_t IT_0219 = IT_0202*IT_0218;
    const ccomplex_t IT_0220 = (-6)*IT_0219;
    const ccomplex_t IT_0221 = 12*s_13;
    const ccomplex_t IT_0222 = IT_0220 + IT_0221;
    const ccomplex_t IT_0223 = pow(m_N_2, 2);
    const ccomplex_t IT_0224 = cpow((-2)*s_12 + IT_0026 + -IT_0027 + -IT_0223 
      + -reg_prop, -1);
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0226 = e_em*IT_0003;
    const ccomplex_t IT_0227 = e_em*IT_0007;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*(IT_0226 + 3*IT_0227);
    const ccomplex_t IT_0229 = (-0.166666666666667)*IT_0228;
    const ccomplex_t IT_0230 = IT_0055*IT_0229;
    const ccomplex_t IT_0231 = IT_0224*IT_0225*IT_0230;
    const ccomplex_t IT_0232 = m_N_2*IT_0039;
    const ccomplex_t IT_0233 = m_N_2*IT_0083;
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = IT_0232 + IT_0234;
    const ccomplex_t IT_0236 = IT_0147*IT_0235;
    const ccomplex_t IT_0237 = m_N_3*IT_0067;
    const ccomplex_t IT_0238 = m_N_2*IT_0156;
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = IT_0237 + IT_0239;
    const ccomplex_t IT_0241 = IT_0170*IT_0240;
    const ccomplex_t IT_0242 = m_N_4*IT_0013;
    const ccomplex_t IT_0243 = m_N_2*IT_0113;
    const ccomplex_t IT_0244 = -IT_0243;
    const ccomplex_t IT_0245 = IT_0242 + IT_0244;
    const ccomplex_t IT_0246 = IT_0127*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0248 = 0.333333333333333*IT_0247;
    const ccomplex_t IT_0249 = (-0.333333333333333)*IT_0052;
    const ccomplex_t IT_0250 = IT_0248*IT_0249;
    const ccomplex_t IT_0251 = IT_0224*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = m_b*IT_0252;
    const ccomplex_t IT_0254 = -IT_0231 + (-0.5)*IT_0236 + -IT_0241 + -IT_0246
       + IT_0253;
    const ccomplex_t IT_0255 = m_N_1*IT_0092;
    const ccomplex_t IT_0256 = m_N_2*IT_0181;
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = IT_0255 + IT_0257;
    const ccomplex_t IT_0259 = IT_0199*IT_0258;
    const ccomplex_t IT_0260 = (-0.5)*IT_0236 + IT_0259;
    const ccomplex_t IT_0261 = IT_0254 + IT_0260;
    const ccomplex_t IT_0262 = conj(IT_0254) + conj(IT_0260);
    const ccomplex_t IT_0263 = IT_0145*IT_0248;
    const ccomplex_t IT_0264 = IT_0224*IT_0225*IT_0263;
    const ccomplex_t IT_0265 = IT_0145*IT_0229;
    const ccomplex_t IT_0266 = IT_0224*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = m_b*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*(0.666666666666667*IT_0042 +
       IT_0045);
    const ccomplex_t IT_0270 = 1.5*IT_0269;
    const ccomplex_t IT_0271 = 0.666666666666667*IT_0270;
    const ccomplex_t IT_0272 = 1.5*IT_0271;
    const ccomplex_t IT_0273 = 0.666666666666667*IT_0272;
    const ccomplex_t IT_0274 = 1.5*IT_0273;
    const ccomplex_t IT_0275 = 0.666666666666667*IT_0274;
    const ccomplex_t IT_0276 = 1.5*IT_0275;
    const ccomplex_t IT_0277 = (-0.333333333333333)*IT_0276;
    const ccomplex_t IT_0278 = IT_0056*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = IT_0235*IT_0279;
    const ccomplex_t IT_0281 = m_N_3*IT_0156;
    const ccomplex_t IT_0282 = m_N_2*IT_0067;
    const ccomplex_t IT_0283 = -IT_0282;
    const ccomplex_t IT_0284 = IT_0281 + IT_0283;
    const ccomplex_t IT_0285 = IT_0080*IT_0284;
    const ccomplex_t IT_0286 = m_N_4*IT_0113;
    const ccomplex_t IT_0287 = m_N_2*IT_0013;
    const ccomplex_t IT_0288 = -IT_0287;
    const ccomplex_t IT_0289 = IT_0286 + IT_0288;
    const ccomplex_t IT_0290 = IT_0030*IT_0289;
    const ccomplex_t IT_0291 = IT_0264 + -IT_0268 + (-0.5)*IT_0280 + IT_0285 +
       IT_0290;
    const ccomplex_t IT_0292 = IT_0058*IT_0235;
    const ccomplex_t IT_0293 = m_N_1*IT_0181;
    const ccomplex_t IT_0294 = m_N_2*IT_0092;
    const ccomplex_t IT_0295 = -IT_0294;
    const ccomplex_t IT_0296 = IT_0293 + IT_0295;
    const ccomplex_t IT_0297 = IT_0103*IT_0296;
    const ccomplex_t IT_0298 = (-0.5)*IT_0292 + -IT_0297;
    const ccomplex_t IT_0299 = IT_0291 + IT_0298;
    const ccomplex_t IT_0300 = conj(IT_0291) + conj(IT_0298);
    const ccomplex_t IT_0301 = m_b*m_N_2;
    const ccomplex_t IT_0302 = pow(s_14, 2);
    const ccomplex_t IT_0303 = IT_0202*IT_0302;
    const ccomplex_t IT_0304 = -IT_0223;
    const ccomplex_t IT_0305 = IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = IT_0301*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0006)*U_sb_00*conj(U_sb_01) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0308 = (-0.666666666666667)*IT_0307;
    const ccomplex_t IT_0309 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0310 = IT_0006*IT_0309;
    const ccomplex_t IT_0311 = 1.4142135623731*IT_0310;
    const ccomplex_t IT_0312 = m_b*conj(N_d2)*e_em*IT_0020*U_sb_11;
    const ccomplex_t IT_0313 = IT_0019*IT_0312;
    const ccomplex_t IT_0314 = 1.4142135623731*IT_0313;
    const ccomplex_t IT_0315 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0316 = IT_0001*IT_0315;
    const ccomplex_t IT_0317 = 1.4142135623731*IT_0316;
    const ccomplex_t IT_0318 = (0 + _Complex_I*1)*(IT_0311 + -IT_0314 + (
      -0.333333333333333)*IT_0317);
    const ccomplex_t IT_0319 = (-3)*IT_0318;
    const ccomplex_t IT_0320 = 0.166666666666667*IT_0319;
    const ccomplex_t IT_0321 = cpow((-2)*s_13 + IT_0026 + IT_0223 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0322 = IT_0320*IT_0321;
    const ccomplex_t IT_0323 = (0 + _Complex_I*1)*IT_0322;
    const ccomplex_t IT_0324 = IT_0308*IT_0323;
    const ccomplex_t IT_0325 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0006)*U_sb_00*conj(U_sb_00) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0326 = (-0.666666666666667)*IT_0325;
    const ccomplex_t IT_0327 = cpow((-2)*s_13 + IT_0026 + IT_0223 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0328 = IT_0145*IT_0327;
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*IT_0328;
    const ccomplex_t IT_0330 = IT_0326*IT_0329;
    const ccomplex_t IT_0331 = (-2)*IT_0128 + -IT_0148 + (-2)*IT_0171 + 
      -IT_0324 + -IT_0330;
    const ccomplex_t IT_0332 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0333 = IT_0001*IT_0332;
    const ccomplex_t IT_0334 = 1.4142135623731*IT_0333;
    const ccomplex_t IT_0335 = m_b*N_d2*e_em*IT_0020*U_sb_01;
    const ccomplex_t IT_0336 = IT_0019*IT_0335;
    const ccomplex_t IT_0337 = 1.4142135623731*IT_0336;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*(IT_0334 + 1.5*IT_0337);
    const ccomplex_t IT_0339 = (-0.333333333333333)*IT_0338;
    const ccomplex_t IT_0340 = IT_0321*IT_0339;
    const ccomplex_t IT_0341 = (0 + _Complex_I*1)*IT_0340;
    const ccomplex_t IT_0342 = IT_0308*IT_0341;
    const ccomplex_t IT_0343 = IT_0056*IT_0249;
    const ccomplex_t IT_0344 = (0 + _Complex_I*1)*IT_0343;
    const ccomplex_t IT_0345 = IT_0039*IT_0344;
    const ccomplex_t IT_0346 = IT_0055*IT_0327;
    const ccomplex_t IT_0347 = (0 + _Complex_I*1)*IT_0346;
    const ccomplex_t IT_0348 = IT_0326*IT_0347;
    const ccomplex_t IT_0349 = 2*IT_0031 + 2*IT_0081 + IT_0342 + IT_0345 +
       IT_0348;
    const ccomplex_t IT_0350 = 6*IT_0349;
    const ccomplex_t IT_0351 = IT_0084 + 2*IT_0104;
    const ccomplex_t IT_0352 = (-6)*IT_0351;
    const ccomplex_t IT_0353 = IT_0350 + IT_0352;
    const ccomplex_t IT_0354 = 6*conj(IT_0349);
    const ccomplex_t IT_0355 = (-6)*conj(IT_0351);
    const ccomplex_t IT_0356 = IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = -IT_0173 + (-2)*IT_0200;
    const ccomplex_t IT_0358 = 6*IT_0351;
    const ccomplex_t IT_0359 = (-6)*IT_0349;
    const ccomplex_t IT_0360 = IT_0358 + IT_0359;
    const ccomplex_t IT_0361 = 6*conj(IT_0351);
    const ccomplex_t IT_0362 = (-6)*conj(IT_0349);
    const ccomplex_t IT_0363 = IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = m_b*IT_0305;
    const ccomplex_t IT_0365 = 6*IT_0331;
    const ccomplex_t IT_0366 = (-6)*IT_0357;
    const ccomplex_t IT_0367 = IT_0365 + IT_0366;
    const ccomplex_t IT_0368 = 6*conj(IT_0331);
    const ccomplex_t IT_0369 = (-6)*conj(IT_0357);
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = s_14*s_34*IT_0202;
    const ccomplex_t IT_0372 = -IT_0371;
    const ccomplex_t IT_0373 = s_13 + IT_0372;
    const ccomplex_t IT_0374 = s_13*IT_0373;
    const ccomplex_t IT_0375 = (-2)*IT_0252 + -IT_0342 + -IT_0348;
    const ccomplex_t IT_0376 = 2*IT_0267 + IT_0324 + IT_0330;
    const ccomplex_t IT_0377 = 6*IT_0357;
    const ccomplex_t IT_0378 = (-6)*IT_0331;
    const ccomplex_t IT_0379 = IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = 6*conj(IT_0357);
    const ccomplex_t IT_0381 = (-6)*conj(IT_0331);
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = s_13*s_14;
    const ccomplex_t IT_0384 = s_34*IT_0223;
    const ccomplex_t IT_0385 = -IT_0384;
    const ccomplex_t IT_0386 = IT_0383 + IT_0385;
    const ccomplex_t IT_0387 = -IT_0252;
    const ccomplex_t IT_0388 = conj(IT_0105) + conj(IT_0387);
    const ccomplex_t IT_0389 = IT_0105 + IT_0387;
    const ccomplex_t IT_0390 = conj(IT_0201) + conj(IT_0267);
    const ccomplex_t IT_0391 = IT_0201 + IT_0267;
    const ccomplex_t IT_0392 = s_34*IT_0209;
    const ccomplex_t IT_0393 = IT_0202*IT_0392;
    const ccomplex_t IT_0394 = (-0.25)*IT_0393;
    const ccomplex_t IT_0395 = s_34 + IT_0394;
    const ccomplex_t IT_0396 = m_N_2*IT_0395;
    const ccomplex_t IT_0397 = (-4)*IT_0396;
    const ccomplex_t IT_0398 = 6*IT_0254;
    const ccomplex_t IT_0399 = 6*IT_0260;
    const ccomplex_t IT_0400 = IT_0398 + IT_0399;
    const ccomplex_t IT_0401 = 6*conj(IT_0254);
    const ccomplex_t IT_0402 = 6*conj(IT_0260);
    const ccomplex_t IT_0403 = IT_0401 + IT_0402;
    const ccomplex_t IT_0404 = 6*IT_0291;
    const ccomplex_t IT_0405 = 6*IT_0298;
    const ccomplex_t IT_0406 = IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = 6*conj(IT_0291);
    const ccomplex_t IT_0408 = 6*conj(IT_0298);
    const ccomplex_t IT_0409 = IT_0407 + IT_0408;
    const ccomplex_t IT_0410 = (-6)*IT_0254;
    const ccomplex_t IT_0411 = (-6)*IT_0260;
    const ccomplex_t IT_0412 = IT_0410 + IT_0411;
    const ccomplex_t IT_0413 = (-6)*conj(IT_0254);
    const ccomplex_t IT_0414 = (-6)*conj(IT_0260);
    const ccomplex_t IT_0415 = IT_0413 + IT_0414;
    const ccomplex_t IT_0416 = (-6)*IT_0291;
    const ccomplex_t IT_0417 = (-6)*IT_0298;
    const ccomplex_t IT_0418 = IT_0416 + IT_0417;
    const ccomplex_t IT_0419 = (-6)*conj(IT_0291);
    const ccomplex_t IT_0420 = (-6)*conj(IT_0298);
    const ccomplex_t IT_0421 = IT_0419 + IT_0420;
    const ccomplex_t IT_0422 = s_14*IT_0209;
    const ccomplex_t IT_0423 = IT_0202*IT_0422;
    const ccomplex_t IT_0424 = 0.5*IT_0423;
    const ccomplex_t IT_0425 = s_14 + IT_0424;
    const ccomplex_t IT_0426 = m_b*IT_0425;
    const ccomplex_t IT_0427 = 2*IT_0426;
    const ccomplex_t IT_0428 = IT_0209*IT_0301;
    const ccomplex_t IT_0429 = 3*IT_0428;
    const ccomplex_t IT_0430 = 6*IT_0387;
    const ccomplex_t IT_0431 = 6*conj(IT_0387);
    const ccomplex_t IT_0432 = (-18)*conj(IT_0254);
    const ccomplex_t IT_0433 = s_13*IT_0305;
    const ccomplex_t IT_0434 = IT_0301*IT_0373;
    const ccomplex_t IT_0435 = 6*IT_0375;
    const ccomplex_t IT_0436 = 6*conj(IT_0375);
    const ccomplex_t IT_0437 = (-6)*IT_0375;
    const ccomplex_t IT_0438 = (-6)*conj(IT_0375);
    const ccomplex_t IT_0439 = pow(s_34, 2);
    const ccomplex_t IT_0440 = IT_0202*IT_0439;
    const ccomplex_t IT_0441 = -IT_0440;
    const ccomplex_t IT_0442 = IT_0026 + IT_0441;
    const ccomplex_t IT_0443 = IT_0301*IT_0442;
    const ccomplex_t IT_0444 = s_13*IT_0442;
    const ccomplex_t IT_0445 = (-6)*IT_0376;
    const ccomplex_t IT_0446 = m_b*IT_0373;
    const ccomplex_t IT_0447 = (-6)*conj(IT_0376);
    const ccomplex_t IT_0448 = s_14*IT_0026;
    const ccomplex_t IT_0449 = s_13*s_34;
    const ccomplex_t IT_0450 = -IT_0449;
    const ccomplex_t IT_0451 = IT_0448 + IT_0450;
    const ccomplex_t IT_0452 = conj(IT_0172)*IT_0376 + IT_0172*conj(IT_0376) +
       0.166666666666667*conj(IT_0082)*IT_0435 + 0.166666666666667*IT_0082
      *IT_0436 + 0.166666666666667*IT_0388*IT_0437 + 0.166666666666667*IT_0389
      *IT_0438 + 0.166666666666667*IT_0390*IT_0445 + 0.166666666666667*IT_0391
      *IT_0447;
    const ccomplex_t IT_0453 = m_N_2*IT_0373;
    const ccomplex_t IT_0454 = m_N_2*IT_0442;
    const ccomplex_t IT_0455 = (-3)*IT_0428;
    const ccomplex_t IT_0456 = 6*IT_0082;
    const ccomplex_t IT_0457 = (-0.166666666666667)*IT_0456;
    const ccomplex_t IT_0458 = IT_0105 + IT_0457;
    const ccomplex_t IT_0459 = 6*conj(IT_0082);
    const ccomplex_t IT_0460 = (-0.166666666666667)*IT_0459;
    const ccomplex_t IT_0461 = conj(IT_0105) + IT_0460;
    const ccomplex_t IT_0462 = (-0.166666666666667)*conj(IT_0267);
    const ccomplex_t IT_0463 = IT_0206 + IT_0214;
    const ccomplex_t IT_0464 = IT_0208 + IT_0213;
    const ccomplex_t IT_0465 = (-0.25)*IT_0423;
    const ccomplex_t IT_0466 = s_14 + IT_0465;
    const ccomplex_t IT_0467 = m_b*IT_0466;
    const ccomplex_t IT_0468 = (-4)*IT_0467;
    const ccomplex_t IT_0469 = (-0.166666666666667)*IT_0468;
    const ccomplex_t IT_0470 = 0.5*IT_0393;
    const ccomplex_t IT_0471 = s_34 + IT_0470;
    const ccomplex_t IT_0472 = m_N_2*IT_0471;
    const ccomplex_t IT_0473 = 2*IT_0472;
    const ccomplex_t IT_0474 = (-0.166666666666667)*IT_0473;
    const ccomplex_t IT_0475 = (conj(IT_0082)*IT_0105 + IT_0082*conj(IT_0105) 
      + conj(IT_0172)*IT_0201 + IT_0172*conj(IT_0201))*IT_0212 + (IT_0082*conj
      (IT_0082) + IT_0105*conj(IT_0105) + IT_0172*conj(IT_0172) + IT_0201*conj
      (IT_0201))*IT_0216 + IT_0222*(IT_0261*IT_0262 + IT_0299*IT_0300) + IT_0306
      *(conj(IT_0331)*IT_0353 + IT_0331*IT_0356 + conj(IT_0357)*IT_0360 +
       IT_0357*IT_0363) + IT_0364*(IT_0262*IT_0353 + IT_0261*IT_0356 + IT_0300
      *IT_0367 + IT_0299*IT_0370) + IT_0374*(IT_0363*IT_0375 + IT_0360*conj
      (IT_0375) + conj(IT_0376)*IT_0379 + IT_0376*IT_0382) + IT_0386*(conj
      (IT_0082)*IT_0353 + IT_0082*IT_0356 + conj(IT_0172)*IT_0367 + IT_0172
      *IT_0370 + IT_0360*IT_0388 + IT_0363*IT_0389 + IT_0379*IT_0390 + IT_0382
      *IT_0391) + IT_0397*(conj(IT_0172)*IT_0400 + IT_0172*IT_0403 + conj
      (IT_0082)*IT_0406 + IT_0082*IT_0409 + conj(IT_0201)*IT_0412 + IT_0201
      *IT_0415 + conj(IT_0105)*IT_0418 + IT_0105*IT_0421) + (conj(IT_0082)
      *IT_0400 + IT_0082*IT_0403 + conj(IT_0172)*IT_0406 + IT_0172*IT_0409 +
       conj(IT_0105)*IT_0412 + IT_0105*IT_0415 + conj(IT_0201)*IT_0418 + IT_0201
      *IT_0421)*IT_0427 + IT_0429*(conj(IT_0172)*IT_0430 + IT_0172*IT_0431) +
       IT_0301*(IT_0261*((-18)*conj(IT_0291) + (-18)*conj(IT_0298)) + IT_0299*((
      -18)*conj(IT_0260) + IT_0432)) + (conj(IT_0349)*IT_0350 + conj(IT_0351)
      *IT_0360 + IT_0351*IT_0362 + conj(IT_0331)*IT_0365 + conj(IT_0357)*IT_0379
       + IT_0357*IT_0381)*IT_0433 + IT_0434*(IT_0363*IT_0376 + IT_0360*conj
      (IT_0376) + conj(IT_0357)*IT_0435 + IT_0357*IT_0436 + conj(IT_0331)
      *IT_0437 + IT_0331*IT_0438) + (conj(IT_0376)*IT_0437 + IT_0376*IT_0438)
      *IT_0443 + IT_0444*(conj(IT_0375)*IT_0437 + conj(IT_0376)*IT_0445) +
       IT_0446*(IT_0262*IT_0437 + IT_0261*IT_0438 + IT_0300*IT_0445 + IT_0299
      *IT_0447) + 6*IT_0451*IT_0452 + (IT_0300*IT_0360 + IT_0299*IT_0363 + conj
      (IT_0260)*IT_0379 + IT_0260*IT_0382 + conj(IT_0357)*IT_0398 + IT_0357
      *IT_0401 + conj(IT_0331)*IT_0410 + IT_0331*IT_0413)*IT_0453 + (conj
      (IT_0376)*IT_0410 + IT_0376*IT_0413 + IT_0300*IT_0437 + IT_0299*IT_0438 +
       conj(IT_0260)*IT_0445 + IT_0260*IT_0447)*IT_0454 + (-6)*IT_0455*((conj
      (IT_0082) + -conj(IT_0105))*IT_0201 + (IT_0082 + -IT_0105)*conj(IT_0201) +
       (-0.166666666666667)*IT_0267*IT_0431 + conj(IT_0172)*IT_0458 + IT_0172
      *IT_0461 + IT_0430*IT_0462) + (-6)*conj(IT_0267)*(IT_0429*IT_0458 + -1./6
      *IT_0201*IT_0463 + -1./6*IT_0172*IT_0464 + IT_0406*IT_0469 + IT_0400
      *IT_0474) + (-6)*IT_0267*(IT_0429*IT_0461 + IT_0216*IT_0462 + -1./6*conj
      (IT_0201)*IT_0463 + -1./6*conj(IT_0172)*IT_0464 + IT_0409*IT_0469 +
       IT_0403*IT_0474) + (-6)*conj(IT_0387)*((-0.166666666666667)*IT_0216
      *IT_0387 + IT_0201*IT_0429 + (-0.166666666666667)*IT_0105*IT_0463 + (
      -0.166666666666667)*IT_0082*IT_0464 + IT_0400*IT_0469 + IT_0406*IT_0474) +
       (-6)*IT_0387*(conj(IT_0201)*IT_0429 + (-0.166666666666667)*conj(IT_0105)
      *IT_0463 + (-0.166666666666667)*conj(IT_0082)*IT_0464 + IT_0403*IT_0469 +
       IT_0409*IT_0474);
    return create_ccomplex_return(IT_0475);
}

