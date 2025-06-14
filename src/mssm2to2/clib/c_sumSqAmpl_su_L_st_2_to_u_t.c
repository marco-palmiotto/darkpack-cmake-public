#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_st_2_to_u_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_st_2_to_u_t(
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
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = IT_0002*IT_0004;
    const ccomplex_t IT_0016 = m_t*conj(N_u1)*e_em*IT_0000*U_st_11;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0011 + 3*IT_0014 + 3
      *IT_0018);
    const ccomplex_t IT_0020 = 0.166666666666667*IT_0019;
    const ccomplex_t IT_0021 = IT_0006*IT_0020;
    const ccomplex_t IT_0022 = pow(m_u, 2);
    const ccomplex_t IT_0023 = pow(m_su_L, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0030 = IT_0008*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0033 = IT_0004*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = m_t*conj(N_u3)*e_em*IT_0000*U_st_11;
    const ccomplex_t IT_0036 = IT_0015*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0031 + 3*IT_0034 + 3
      *IT_0037);
    const ccomplex_t IT_0039 = 0.166666666666667*IT_0038;
    const ccomplex_t IT_0040 = IT_0028*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0047 = IT_0008*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0050 = IT_0004*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = m_t*conj(N_u4)*e_em*IT_0000*U_st_11;
    const ccomplex_t IT_0053 = IT_0015*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0048 + 3*IT_0051 + 3
      *IT_0054);
    const ccomplex_t IT_0056 = 0.166666666666667*IT_0055;
    const ccomplex_t IT_0057 = IT_0045*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0062 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0063 = IT_0008*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0066 = IT_0004*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = m_t*conj(N_u2)*e_em*IT_0000*U_st_11;
    const ccomplex_t IT_0069 = IT_0015*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0064 + 3*IT_0067 + 3
      *IT_0070);
    const ccomplex_t IT_0072 = 0.166666666666667*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = IT_0061*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = -IT_0026 + -IT_0043 + -IT_0060 + -IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = pow(m_sG, 2);
    const ccomplex_t IT_0083 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -IT_0082 +
       reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B1)*e_em;
    const ccomplex_t IT_0087 = IT_0008*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W1)*e_em;
    const ccomplex_t IT_0090 = IT_0004*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + 3*IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0095 = IT_0008*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = m_t*N_u1*e_em*IT_0000*U_st_01;
    const ccomplex_t IT_0098 = IT_0015*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + (-0.75)*IT_0099);
    const ccomplex_t IT_0101 = 0.666666666666667*IT_0100;
    const ccomplex_t IT_0102 = IT_0093*IT_0101;
    const ccomplex_t IT_0103 = IT_0024*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0106 = IT_0008*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = m_t*N_u2*e_em*IT_0000*U_st_01;
    const ccomplex_t IT_0109 = IT_0015*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0107 + (-0.75)*IT_0110);
    const ccomplex_t IT_0112 = 0.666666666666667*IT_0111;
    const ccomplex_t IT_0113 = conj(N_B2)*e_em;
    const ccomplex_t IT_0114 = IT_0008*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = conj(N_W2)*e_em;
    const ccomplex_t IT_0117 = IT_0004*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + 3*IT_0118);
    const ccomplex_t IT_0120 = 0.166666666666667*IT_0119;
    const ccomplex_t IT_0121 = IT_0112*IT_0120;
    const ccomplex_t IT_0122 = IT_0061*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = conj(N_B3)*e_em;
    const ccomplex_t IT_0125 = IT_0008*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = conj(N_W3)*e_em;
    const ccomplex_t IT_0128 = IT_0004*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + 3*IT_0129);
    const ccomplex_t IT_0131 = 0.166666666666667*IT_0130;
    const ccomplex_t IT_0132 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0133 = IT_0008*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = m_t*N_u3*e_em*IT_0000*U_st_01;
    const ccomplex_t IT_0136 = IT_0015*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0134 + (-0.75)*IT_0137);
    const ccomplex_t IT_0139 = 0.666666666666667*IT_0138;
    const ccomplex_t IT_0140 = IT_0131*IT_0139;
    const ccomplex_t IT_0141 = IT_0041*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = conj(N_B4)*e_em;
    const ccomplex_t IT_0144 = IT_0008*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = conj(N_W4)*e_em;
    const ccomplex_t IT_0147 = IT_0004*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + 3*IT_0148);
    const ccomplex_t IT_0150 = 0.166666666666667*IT_0149;
    const ccomplex_t IT_0151 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0152 = IT_0008*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = m_t*N_u4*e_em*IT_0000*U_st_01;
    const ccomplex_t IT_0155 = IT_0015*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0153 + (-0.75)*IT_0156);
    const ccomplex_t IT_0158 = 0.666666666666667*IT_0157;
    const ccomplex_t IT_0159 = IT_0150*IT_0158;
    const ccomplex_t IT_0160 = IT_0058*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = 0.166666666666667*IT_0085 + -IT_0104 + -IT_0123
       + -IT_0142 + -IT_0161;
    const ccomplex_t IT_0163 = (-0.5)*IT_0085;
    const ccomplex_t IT_0164 = m_t*m_u*IT_0023;
    const ccomplex_t IT_0165 = s_34*IT_0023;
    const ccomplex_t IT_0166 = (-18)*IT_0165;
    const ccomplex_t IT_0167 = s_13*s_14;
    const ccomplex_t IT_0168 = 36*IT_0167;
    const ccomplex_t IT_0169 = IT_0166 + IT_0168;
    const ccomplex_t IT_0170 = m_t*m_u;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0172 = IT_0079*IT_0171;
    const ccomplex_t IT_0173 = IT_0083*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = m_sG*IT_0174;
    const ccomplex_t IT_0176 = (-0.5)*IT_0175;
    const ccomplex_t IT_0177 = m_u*IT_0085;
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = IT_0006*IT_0101;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0181 = IT_0024*IT_0179*IT_0180;
    const ccomplex_t IT_0182 = IT_0045*IT_0158;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0184 = IT_0058*IT_0182*IT_0183;
    const ccomplex_t IT_0185 = m_u*IT_0104;
    const ccomplex_t IT_0186 = m_u*IT_0123;
    const ccomplex_t IT_0187 = m_u*IT_0142;
    const ccomplex_t IT_0188 = m_u*IT_0161;
    const ccomplex_t IT_0189 = 0.166666666666667*IT_0177 + IT_0181 + IT_0184 +
       -IT_0185 + -IT_0186 + -IT_0187 + -IT_0188;
    const ccomplex_t IT_0190 = IT_0074*IT_0112;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0192 = IT_0061*IT_0190*IT_0191;
    const ccomplex_t IT_0193 = IT_0028*IT_0139;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0195 = IT_0041*IT_0193*IT_0194;
    const ccomplex_t IT_0196 = -IT_0192 + -IT_0195;
    const ccomplex_t IT_0197 = IT_0020*IT_0093;
    const ccomplex_t IT_0198 = IT_0024*IT_0180*IT_0197;
    const ccomplex_t IT_0199 = IT_0056*IT_0150;
    const ccomplex_t IT_0200 = IT_0058*IT_0183*IT_0199;
    const ccomplex_t IT_0201 = IT_0072*IT_0120;
    const ccomplex_t IT_0202 = IT_0061*IT_0191*IT_0201;
    const ccomplex_t IT_0203 = 0.166666666666667*IT_0175 + -IT_0198 + -IT_0200
       + -IT_0202;
    const ccomplex_t IT_0204 = 6*IT_0178;
    const ccomplex_t IT_0205 = 18*IT_0189;
    const ccomplex_t IT_0206 = (-18)*IT_0196;
    const ccomplex_t IT_0207 = IT_0204 + IT_0205 + IT_0206;
    const ccomplex_t IT_0208 = 6*conj(IT_0178);
    const ccomplex_t IT_0209 = 18*conj(IT_0189);
    const ccomplex_t IT_0210 = (-18)*conj(IT_0196);
    const ccomplex_t IT_0211 = IT_0208 + IT_0209 + IT_0210;
    const ccomplex_t IT_0212 = m_u*IT_0026;
    const ccomplex_t IT_0213 = m_u*IT_0043;
    const ccomplex_t IT_0214 = m_u*IT_0060;
    const ccomplex_t IT_0215 = IT_0039*IT_0131;
    const ccomplex_t IT_0216 = IT_0041*IT_0194*IT_0215;
    const ccomplex_t IT_0217 = m_u*IT_0077;
    const ccomplex_t IT_0218 = -IT_0212 + -IT_0213 + -IT_0214 + IT_0216 + 
      -IT_0217;
    const ccomplex_t IT_0219 = 18*IT_0196;
    const ccomplex_t IT_0220 = (-18)*IT_0189;
    const ccomplex_t IT_0221 = (-6)*IT_0178;
    const ccomplex_t IT_0222 = IT_0219 + IT_0220 + IT_0221;
    const ccomplex_t IT_0223 = 18*conj(IT_0196);
    const ccomplex_t IT_0224 = (-18)*conj(IT_0189);
    const ccomplex_t IT_0225 = (-6)*conj(IT_0178);
    const ccomplex_t IT_0226 = IT_0223 + IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = conj(IT_0176)*(18*IT_0178 + 6*IT_0189 + (-6)
      *IT_0196) + IT_0176*(18*conj(IT_0178) + 6*conj(IT_0189) + (-6)*conj
      (IT_0196)) + conj(IT_0203)*IT_0207 + IT_0203*IT_0211 + conj(IT_0218)
      *IT_0222 + IT_0218*IT_0226;
    const ccomplex_t IT_0228 = 6*s_34;
    const ccomplex_t IT_0229 = 18*s_34;
    const ccomplex_t IT_0230 = (-18)*s_34;
    const ccomplex_t IT_0231 = (-6)*s_34;
    const ccomplex_t IT_0232 = 12*IT_0167;
    const ccomplex_t IT_0233 = (-6)*IT_0165;
    const ccomplex_t IT_0234 = m_u*s_14;
    const ccomplex_t IT_0235 = (-3)*IT_0218;
    const ccomplex_t IT_0236 = 3*conj(IT_0203);
    const ccomplex_t IT_0237 = (-3)*conj(IT_0218);
    const ccomplex_t IT_0238 = conj(IT_0163)*(IT_0178 + 0.333333333333333
      *IT_0189 + (-0.333333333333333)*IT_0196) + IT_0163*(conj(IT_0178) +
       0.333333333333333*conj(IT_0189) + (-0.333333333333333)*conj(IT_0196)) + (
      -0.0555555555555556)*conj(IT_0162)*IT_0222 + (-0.0555555555555556)*IT_0162
      *IT_0226 + (-0.333333333333333)*conj(IT_0078)*(IT_0176 + 3*IT_0203 +
       IT_0235) + (-0.333333333333333)*IT_0078*(conj(IT_0176) + IT_0236 +
       IT_0237);
    const ccomplex_t IT_0239 = m_t*s_13;
    const ccomplex_t IT_0240 = (-0.0555555555555556)*conj(IT_0078)*IT_0207 + (
      -0.0555555555555556)*IT_0078*IT_0211 + conj(IT_0163)*(IT_0176 +
       0.333333333333333*IT_0203 + (-0.333333333333333)*IT_0218) + IT_0163*(conj
      (IT_0176) + 0.333333333333333*conj(IT_0203) + (-0.333333333333333)*conj
      (IT_0218)) + 0.333333333333333*conj(IT_0162)*(IT_0176 + 3*IT_0203 +
       IT_0235) + 0.333333333333333*IT_0162*(conj(IT_0176) + IT_0236 + IT_0237);
    const ccomplex_t IT_0241 = (conj(IT_0078)*((-18)*IT_0162 + (-6)*IT_0163) +
       IT_0078*((-18)*conj(IT_0162) + (-6)*conj(IT_0163)))*IT_0164 + (IT_0078
      *conj(IT_0078) + IT_0162*conj(IT_0162) + IT_0163*conj(IT_0163))*IT_0169 +
       IT_0170*IT_0227 + (conj(IT_0178)*IT_0189 + IT_0178*conj(IT_0189) + conj
      (IT_0176)*IT_0203 + IT_0176*conj(IT_0203))*IT_0228 + (IT_0176*conj(IT_0176
      ) + IT_0178*conj(IT_0178) + IT_0189*conj(IT_0189) + IT_0196*conj(IT_0196) 
      + IT_0203*conj(IT_0203) + IT_0218*conj(IT_0218))*IT_0229 + (conj(IT_0189)
      *IT_0196 + IT_0189*conj(IT_0196) + conj(IT_0203)*IT_0218 + IT_0203*conj
      (IT_0218))*IT_0230 + (conj(IT_0178)*IT_0196 + IT_0178*conj(IT_0196) + conj
      (IT_0176)*IT_0218 + IT_0176*conj(IT_0218))*IT_0231 + (conj(IT_0162)
      *IT_0163 + IT_0162*conj(IT_0163))*(IT_0232 + IT_0233) + (-18)*IT_0234
      *IT_0238 + (-18)*IT_0239*IT_0240;
    return create_ccomplex_return(IT_0241);
}

