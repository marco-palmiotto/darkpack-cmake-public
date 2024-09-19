#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_st_1_to_u_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_st_1_to_u_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
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
    const creal_t m_su_L = param->m_su_L;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_su_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = IT_0005*IT_0007;
    const ccomplex_t IT_0014 = m_t*conj(N_u1)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = cos(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0012 + IT_0016 +
       0.333333333333333*IT_0021);
    const ccomplex_t IT_0023 = 3*IT_0022;
    const ccomplex_t IT_0024 = 0.166666666666667*IT_0023;
    const ccomplex_t IT_0025 = IT_0009*IT_0024;
    const ccomplex_t IT_0026 = IT_0002*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0032 = IT_0018*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0035 = IT_0007*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = m_t*conj(N_u2)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0038 = IT_0013*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0033 + 3*IT_0036 + 3
      *IT_0039);
    const ccomplex_t IT_0041 = 0.166666666666667*IT_0040;
    const ccomplex_t IT_0042 = IT_0030*IT_0041;
    const ccomplex_t IT_0043 = IT_0028*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0049 = IT_0007*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = m_t*conj(N_u3)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0052 = IT_0013*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0055 = IT_0018*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0050 + IT_0053 +
       0.333333333333333*IT_0056);
    const ccomplex_t IT_0058 = 3*IT_0057;
    const ccomplex_t IT_0059 = 0.166666666666667*IT_0058;
    const ccomplex_t IT_0060 = IT_0047*IT_0059;
    const ccomplex_t IT_0061 = IT_0045*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0067 = IT_0007*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = m_t*conj(N_u4)*e_em*IT_0003*U_st_10;
    const ccomplex_t IT_0070 = IT_0013*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0073 = IT_0018*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0068 + IT_0071 +
       0.333333333333333*IT_0074);
    const ccomplex_t IT_0076 = 3*IT_0075;
    const ccomplex_t IT_0077 = 0.166666666666667*IT_0076;
    const ccomplex_t IT_0078 = IT_0065*IT_0077;
    const ccomplex_t IT_0079 = IT_0063*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = -IT_0027 + -IT_0044 + -IT_0062 + -IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = pow(m_sG, 2);
    const ccomplex_t IT_0086 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0085 +
       reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = conj(N_B1)*e_em;
    const ccomplex_t IT_0091 = IT_0018*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W1)*e_em;
    const ccomplex_t IT_0094 = IT_0007*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0092 + 3*IT_0095);
    const ccomplex_t IT_0097 = 0.166666666666667*IT_0096;
    const ccomplex_t IT_0098 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0099 = IT_0018*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = m_t*N_u1*e_em*IT_0003*U_st_00;
    const ccomplex_t IT_0102 = IT_0013*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + (-0.75)*IT_0103);
    const ccomplex_t IT_0105 = 0.666666666666667*IT_0104;
    const ccomplex_t IT_0106 = IT_0097*IT_0105;
    const ccomplex_t IT_0107 = IT_0002*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B2)*e_em;
    const ccomplex_t IT_0110 = IT_0018*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_W2)*e_em;
    const ccomplex_t IT_0113 = IT_0007*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + 3*IT_0114);
    const ccomplex_t IT_0116 = 0.166666666666667*IT_0115;
    const ccomplex_t IT_0117 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0118 = IT_0018*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = m_t*N_u2*e_em*IT_0003*U_st_00;
    const ccomplex_t IT_0121 = IT_0013*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0119 + (-0.75)*IT_0122);
    const ccomplex_t IT_0124 = 0.666666666666667*IT_0123;
    const ccomplex_t IT_0125 = IT_0116*IT_0124;
    const ccomplex_t IT_0126 = IT_0028*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = conj(N_B3)*e_em;
    const ccomplex_t IT_0129 = IT_0018*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = conj(N_W3)*e_em;
    const ccomplex_t IT_0132 = IT_0007*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0130 + 3*IT_0133);
    const ccomplex_t IT_0135 = 0.166666666666667*IT_0134;
    const ccomplex_t IT_0136 = m_t*N_u3*e_em*IT_0003*U_st_00;
    const ccomplex_t IT_0137 = IT_0013*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0140 = IT_0018*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0138 + (
      -1.33333333333333)*IT_0141);
    const ccomplex_t IT_0143 = (-0.75)*IT_0142;
    const ccomplex_t IT_0144 = (-1.33333333333333)*IT_0143;
    const ccomplex_t IT_0145 = (-0.75)*IT_0144;
    const ccomplex_t IT_0146 = (-1.33333333333333)*IT_0145;
    const ccomplex_t IT_0147 = (-0.75)*IT_0146;
    const ccomplex_t IT_0148 = (-1.33333333333333)*IT_0147;
    const ccomplex_t IT_0149 = (-0.75)*IT_0148;
    const ccomplex_t IT_0150 = (-1.33333333333333)*IT_0149;
    const ccomplex_t IT_0151 = (-0.75)*IT_0150;
    const ccomplex_t IT_0152 = (-1.33333333333333)*IT_0151;
    const ccomplex_t IT_0153 = (-0.75)*IT_0152;
    const ccomplex_t IT_0154 = (-1.33333333333333)*IT_0153;
    const ccomplex_t IT_0155 = (-0.75)*IT_0154;
    const ccomplex_t IT_0156 = 0.666666666666667*IT_0155;
    const ccomplex_t IT_0157 = IT_0135*IT_0156;
    const ccomplex_t IT_0158 = IT_0045*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = conj(N_B4)*e_em;
    const ccomplex_t IT_0161 = IT_0018*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = conj(N_W4)*e_em;
    const ccomplex_t IT_0164 = IT_0007*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*(IT_0162 + 3*IT_0165);
    const ccomplex_t IT_0167 = 0.166666666666667*IT_0166;
    const ccomplex_t IT_0168 = m_t*N_u4*e_em*IT_0003*U_st_00;
    const ccomplex_t IT_0169 = IT_0013*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0172 = IT_0018*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0170 + (
      -1.33333333333333)*IT_0173);
    const ccomplex_t IT_0175 = (-0.75)*IT_0174;
    const ccomplex_t IT_0176 = 0.666666666666667*IT_0175;
    const ccomplex_t IT_0177 = IT_0167*IT_0176;
    const ccomplex_t IT_0178 = IT_0063*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = 0.166666666666667*IT_0088 + -IT_0108 + -IT_0127
       + -IT_0159 + -IT_0179;
    const ccomplex_t IT_0181 = m_t*m_u*IT_0001;
    const ccomplex_t IT_0182 = s_34*IT_0001;
    const ccomplex_t IT_0183 = (-18)*IT_0182;
    const ccomplex_t IT_0184 = s_13*s_14;
    const ccomplex_t IT_0185 = 36*IT_0184;
    const ccomplex_t IT_0186 = IT_0183 + IT_0185;
    const ccomplex_t IT_0187 = m_t*m_u;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0189 = IT_0082*IT_0188;
    const ccomplex_t IT_0190 = IT_0086*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = m_sG*IT_0191;
    const ccomplex_t IT_0193 = 0.5*IT_0192;
    const ccomplex_t IT_0194 = IT_0030*IT_0124;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0196 = IT_0028*IT_0194*IT_0195;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = 6*IT_0197;
    const ccomplex_t IT_0199 = m_u*IT_0088;
    const ccomplex_t IT_0200 = (-0.5)*IT_0199;
    const ccomplex_t IT_0201 = (-18)*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0203 = IT_0065*IT_0176;
    const ccomplex_t IT_0204 = IT_0063*IT_0202*IT_0203;
    const ccomplex_t IT_0205 = m_u*IT_0127;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0207 = IT_0009*IT_0105;
    const ccomplex_t IT_0208 = IT_0002*IT_0206*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0210 = IT_0047*IT_0156;
    const ccomplex_t IT_0211 = IT_0045*IT_0209*IT_0210;
    const ccomplex_t IT_0212 = m_u*IT_0108;
    const ccomplex_t IT_0213 = m_u*IT_0159;
    const ccomplex_t IT_0214 = m_u*IT_0179;
    const ccomplex_t IT_0215 = 0.166666666666667*IT_0199 + IT_0204 + -IT_0205 
      + IT_0208 + IT_0211 + -IT_0212 + -IT_0213 + -IT_0214;
    const ccomplex_t IT_0216 = (-6)*IT_0215;
    const ccomplex_t IT_0217 = IT_0198 + IT_0201 + IT_0216;
    const ccomplex_t IT_0218 = 6*conj(IT_0197);
    const ccomplex_t IT_0219 = (-18)*conj(IT_0200);
    const ccomplex_t IT_0220 = (-6)*conj(IT_0215);
    const ccomplex_t IT_0221 = IT_0218 + IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = IT_0077*IT_0167;
    const ccomplex_t IT_0223 = IT_0063*IT_0202*IT_0222;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = IT_0024*IT_0097;
    const ccomplex_t IT_0226 = IT_0002*IT_0206*IT_0225;
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = conj(IT_0224) + conj(IT_0227);
    const ccomplex_t IT_0229 = 6*IT_0200;
    const ccomplex_t IT_0230 = 18*IT_0215;
    const ccomplex_t IT_0231 = (-18)*IT_0197;
    const ccomplex_t IT_0232 = IT_0229 + IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = IT_0224 + IT_0227;
    const ccomplex_t IT_0234 = 6*conj(IT_0200);
    const ccomplex_t IT_0235 = 18*conj(IT_0215);
    const ccomplex_t IT_0236 = (-18)*conj(IT_0197);
    const ccomplex_t IT_0237 = IT_0234 + IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = IT_0059*IT_0135;
    const ccomplex_t IT_0239 = IT_0045*IT_0209*IT_0238;
    const ccomplex_t IT_0240 = m_u*IT_0080;
    const ccomplex_t IT_0241 = IT_0041*IT_0116;
    const ccomplex_t IT_0242 = IT_0028*IT_0195*IT_0241;
    const ccomplex_t IT_0243 = m_u*IT_0027;
    const ccomplex_t IT_0244 = m_u*IT_0044;
    const ccomplex_t IT_0245 = m_u*IT_0062;
    const ccomplex_t IT_0246 = (-0.166666666666667)*IT_0192 + IT_0239 + 
      -IT_0240 + IT_0242 + -IT_0243 + -IT_0244 + -IT_0245;
    const ccomplex_t IT_0247 = 18*IT_0197;
    const ccomplex_t IT_0248 = (-18)*IT_0215;
    const ccomplex_t IT_0249 = (-6)*IT_0200;
    const ccomplex_t IT_0250 = IT_0247 + IT_0248 + IT_0249;
    const ccomplex_t IT_0251 = 18*conj(IT_0197);
    const ccomplex_t IT_0252 = (-18)*conj(IT_0215);
    const ccomplex_t IT_0253 = (-6)*conj(IT_0200);
    const ccomplex_t IT_0254 = IT_0251 + IT_0252 + IT_0253;
    const ccomplex_t IT_0255 = m_u*s_14;
    const ccomplex_t IT_0256 = conj(IT_0089)*IT_0217 + IT_0089*IT_0221 + (-6)
      *conj(IT_0081)*(IT_0193 + (-3)*IT_0224 + (-3)*IT_0227 + 3*IT_0246) + (-6)
      *IT_0081*(conj(IT_0193) + (-3)*conj(IT_0224) + (-3)*conj(IT_0227) + 3*conj
      (IT_0246)) + conj(IT_0180)*IT_0250 + IT_0180*IT_0254;
    const ccomplex_t IT_0257 = 6*s_34;
    const ccomplex_t IT_0258 = 18*s_34;
    const ccomplex_t IT_0259 = (-18)*s_34;
    const ccomplex_t IT_0260 = IT_0224*IT_0258;
    const ccomplex_t IT_0261 = IT_0246*IT_0259;
    const ccomplex_t IT_0262 = conj(IT_0246)*IT_0259;
    const ccomplex_t IT_0263 = (-6)*s_34;
    const ccomplex_t IT_0264 = 12*IT_0184;
    const ccomplex_t IT_0265 = (-6)*IT_0182;
    const ccomplex_t IT_0266 = m_t*s_13;
    const ccomplex_t IT_0267 = -IT_0224 + IT_0246;
    const ccomplex_t IT_0268 = -IT_0227;
    const ccomplex_t IT_0269 = -conj(IT_0224);
    const ccomplex_t IT_0270 = conj(IT_0246) + IT_0269;
    const ccomplex_t IT_0271 = -conj(IT_0227);
    const ccomplex_t IT_0272 = conj(IT_0081)*IT_0232 + IT_0081*IT_0237 + 6
      *conj(IT_0180)*(IT_0193 + 3*IT_0267 + 3*IT_0268) + 18*conj(IT_0089)*
      (IT_0193 + 0.333333333333333*IT_0267 + 0.333333333333333*IT_0268) + 6
      *IT_0180*(conj(IT_0193) + 3*IT_0270 + 3*IT_0271) + 18*IT_0089*(conj
      (IT_0193) + 0.333333333333333*IT_0270 + 0.333333333333333*IT_0271);
    const ccomplex_t IT_0273 = (-6)*(conj(IT_0081)*(IT_0089 + 3*IT_0180) +
       IT_0081*(conj(IT_0089) + 3*conj(IT_0180)))*IT_0181 + (IT_0081*conj
      (IT_0081) + IT_0089*conj(IT_0089) + IT_0180*conj(IT_0180))*IT_0186 +
       IT_0187*(conj(IT_0193)*IT_0217 + IT_0193*IT_0221 + IT_0228*IT_0232 +
       IT_0233*IT_0237 + conj(IT_0246)*IT_0250 + IT_0246*IT_0254) + IT_0255
      *IT_0256 + (conj(IT_0200)*IT_0215 + IT_0200*conj(IT_0215) + conj(IT_0193)
      *IT_0246 + IT_0193*conj(IT_0246))*IT_0257 + (IT_0193*conj(IT_0193) +
       IT_0197*conj(IT_0197) + IT_0200*conj(IT_0200) + IT_0215*conj(IT_0215) +
       IT_0227*IT_0228 + IT_0246*conj(IT_0246))*IT_0258 + (conj(IT_0197)*IT_0215
       + IT_0197*conj(IT_0215))*IT_0259 + IT_0228*(IT_0260 + IT_0261) + IT_0233
      *IT_0262 + (conj(IT_0197)*IT_0200 + IT_0197*conj(IT_0200) + IT_0193
      *IT_0228 + conj(IT_0193)*IT_0233)*IT_0263 + (conj(IT_0089)*IT_0180 +
       IT_0089*conj(IT_0180))*(IT_0264 + IT_0265) + IT_0266*IT_0272;
    return create_ccomplex_return(IT_0273);
}

