#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_L_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_L_to_anti_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_sG = param->m_sG;
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
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = pow(m_s, 2);
    const ccomplex_t IT_0009 = cos(alpha);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0003*IT_0005*IT_0008 + (-0.5)*m_W*(IT_0001*IT_0003 + IT_0009*IT_0010)*
      (IT_0005 + 1./3*IT_0004*IT_0012));
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = pow(m_ss_L, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0005*IT_0009;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0005*IT_0008*IT_0009 + (-0.5)*m_W*(IT_0001*IT_0009 + -IT_0003*IT_0010)
      *(IT_0005 + 1./3*IT_0004*IT_0012));
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + (-2)*IT_0015 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0018 + -IT_0026;
    const ccomplex_t IT_0028 = cpow(IT_0011, -1);
    const ccomplex_t IT_0029 = IT_0004*IT_0028;
    const ccomplex_t IT_0030 = e_em*IT_0029;
    const ccomplex_t IT_0031 = IT_0005*IT_0011;
    const ccomplex_t IT_0032 = e_em*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0030 + 3*IT_0032);
    const ccomplex_t IT_0034 = (-0.166666666666667)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*3)*e_em*(IT_0005*IT_0011 + 1./3
      *IT_0004*IT_0028);
    const ccomplex_t IT_0039 = 0.166666666666667*IT_0038;
    const ccomplex_t IT_0040 = m_s*IT_0039;
    const ccomplex_t IT_0041 = 0.333333333333333*IT_0038;
    const ccomplex_t IT_0042 = m_s*IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0040 + IT_0043;
    const ccomplex_t IT_0045 = IT_0037*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0028;
    const ccomplex_t IT_0047 = 0.333333333333333*IT_0046;
    const ccomplex_t IT_0048 = IT_0035*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0044*IT_0049;
    const ccomplex_t IT_0051 = pow(m_sG, 2);
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0008 + IT_0015 + -IT_0051 +
       reg_prop, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0054 = cpow(IT_0053, 2);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0052*IT_0055;
    const ccomplex_t IT_0057 = m_s*IT_0056;
    const ccomplex_t IT_0058 = -IT_0045 + IT_0050 + 0.5*IT_0057;
    const ccomplex_t IT_0059 = (IT_0040 + -IT_0042)*(IT_0037 + -IT_0049);
    const ccomplex_t IT_0060 = (-18)*IT_0008;
    const ccomplex_t IT_0061 = 18*s_34;
    const ccomplex_t IT_0062 = IT_0060 + IT_0061;
    const ccomplex_t IT_0063 = IT_0037*IT_0041;
    const ccomplex_t IT_0064 = (-0.5)*IT_0056 + -IT_0063;
    const ccomplex_t IT_0065 = conj(N_B1)*e_em;
    const ccomplex_t IT_0066 = IT_0028*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W1)*e_em;
    const ccomplex_t IT_0069 = IT_0005*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0067 + (-3)*IT_0070);
    const ccomplex_t IT_0072 = 0.166666666666667*IT_0071;
    const ccomplex_t IT_0073 = N_B1*e_em;
    const ccomplex_t IT_0074 = IT_0028*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = N_W1*e_em;
    const ccomplex_t IT_0077 = IT_0005*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0075 + (-3)*IT_0078);
    const ccomplex_t IT_0080 = 0.166666666666667*IT_0079;
    const ccomplex_t IT_0081 = IT_0072*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_23 + IT_0008 + IT_0015 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = conj(N_B2)*e_em;
    const ccomplex_t IT_0086 = IT_0028*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W2)*e_em;
    const ccomplex_t IT_0089 = IT_0005*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + (-3)*IT_0090);
    const ccomplex_t IT_0092 = 0.166666666666667*IT_0091;
    const ccomplex_t IT_0093 = N_B2*e_em;
    const ccomplex_t IT_0094 = IT_0028*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = N_W2*e_em;
    const ccomplex_t IT_0097 = IT_0005*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0095 + (-3)*IT_0098);
    const ccomplex_t IT_0100 = 0.166666666666667*IT_0099;
    const ccomplex_t IT_0101 = IT_0092*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_23 + IT_0008 + IT_0015 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_23 + IT_0008 + IT_0015 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0106 = N_B3*e_em;
    const ccomplex_t IT_0107 = IT_0028*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = N_W3*e_em;
    const ccomplex_t IT_0110 = IT_0005*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + (-3)*IT_0111);
    const ccomplex_t IT_0113 = 0.166666666666667*IT_0112;
    const ccomplex_t IT_0114 = conj(N_B3)*e_em;
    const ccomplex_t IT_0115 = IT_0028*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = conj(N_W3)*e_em;
    const ccomplex_t IT_0118 = IT_0005*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0116 + (-3)*IT_0119);
    const ccomplex_t IT_0121 = 0.166666666666667*IT_0120;
    const ccomplex_t IT_0122 = IT_0113*IT_0121;
    const ccomplex_t IT_0123 = IT_0105*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = conj(N_B4)*e_em;
    const ccomplex_t IT_0126 = IT_0028*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = conj(N_W4)*e_em;
    const ccomplex_t IT_0129 = IT_0005*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + (-3)*IT_0130);
    const ccomplex_t IT_0132 = 0.166666666666667*IT_0131;
    const ccomplex_t IT_0133 = N_B4*e_em;
    const ccomplex_t IT_0134 = IT_0028*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = N_W4*e_em;
    const ccomplex_t IT_0137 = IT_0005*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0135 + (-3)*IT_0138);
    const ccomplex_t IT_0140 = 0.166666666666667*IT_0139;
    const ccomplex_t IT_0141 = IT_0132*IT_0140;
    const ccomplex_t IT_0142 = cpow((-2)*s_23 + IT_0008 + IT_0015 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0143 = IT_0141*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = 0.166666666666667*IT_0056 + -IT_0084 + -IT_0104
       + -IT_0124 + -IT_0144;
    const ccomplex_t IT_0146 = IT_0041*IT_0049;
    const ccomplex_t IT_0147 = -IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = IT_0149*IT_0151;
    const ccomplex_t IT_0153 = IT_0082*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0158 = (-0.5)*IT_0157;
    const ccomplex_t IT_0159 = IT_0156*IT_0158;
    const ccomplex_t IT_0160 = IT_0142*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0163 = (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0165 = (-0.5)*IT_0164;
    const ccomplex_t IT_0166 = IT_0163*IT_0165;
    const ccomplex_t IT_0167 = IT_0102*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0172 = (-0.5)*IT_0171;
    const ccomplex_t IT_0173 = IT_0170*IT_0172;
    const ccomplex_t IT_0174 = IT_0105*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = -IT_0154 + -IT_0161 + -IT_0168 + -IT_0175;
    const ccomplex_t IT_0177 = s_34*IT_0015;
    const ccomplex_t IT_0178 = (-6)*IT_0177;
    const ccomplex_t IT_0179 = s_13*s_14;
    const ccomplex_t IT_0180 = 12*IT_0179;
    const ccomplex_t IT_0181 = IT_0178 + IT_0180;
    const ccomplex_t IT_0182 = (-18)*IT_0177;
    const ccomplex_t IT_0183 = 36*IT_0179;
    const ccomplex_t IT_0184 = IT_0182 + IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0186 = (-2)*IT_0185;
    const ccomplex_t IT_0187 = cpow(s_12 + IT_0015 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0188 = IT_0185*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = IT_0186*IT_0189;
    const ccomplex_t IT_0191 = 0.25*IT_0190;
    const ccomplex_t IT_0192 = cpow(IT_0191, 2);
    const ccomplex_t IT_0193 = (-36)*IT_0177;
    const ccomplex_t IT_0194 = IT_0008*IT_0015;
    const ccomplex_t IT_0195 = (-36)*IT_0194;
    const ccomplex_t IT_0196 = 72*IT_0179;
    const ccomplex_t IT_0197 = IT_0193 + IT_0195 + IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0199 = 0.666666666666667*IT_0198;
    const ccomplex_t IT_0200 = (-0.333333333333333)*IT_0198;
    const ccomplex_t IT_0201 = IT_0199*IT_0200;
    const ccomplex_t IT_0202 = IT_0187*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = (-0.0833333333333333)*IT_0190 + 0.5*IT_0203;
    const ccomplex_t IT_0205 = (-18)*IT_0058;
    const ccomplex_t IT_0206 = IT_0080*IT_0151;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0208 = IT_0082*IT_0206*IT_0207;
    const ccomplex_t IT_0209 = IT_0100*IT_0165;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0211 = IT_0102*IT_0209*IT_0210;
    const ccomplex_t IT_0212 = IT_0113*IT_0170;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0214 = IT_0105*IT_0212*IT_0213;
    const ccomplex_t IT_0215 = IT_0140*IT_0158;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0217 = IT_0142*IT_0215*IT_0216;
    const ccomplex_t IT_0218 = m_s*IT_0084;
    const ccomplex_t IT_0219 = m_s*IT_0104;
    const ccomplex_t IT_0220 = m_s*IT_0124;
    const ccomplex_t IT_0221 = m_s*IT_0144;
    const ccomplex_t IT_0222 = (-0.166666666666667)*IT_0057 + -IT_0208 + 
      -IT_0211 + -IT_0214 + -IT_0217 + IT_0218 + IT_0219 + IT_0220 + IT_0221;
    const ccomplex_t IT_0223 = (-6)*IT_0222;
    const ccomplex_t IT_0224 = IT_0205 + IT_0223;
    const ccomplex_t IT_0225 = (-18)*conj(IT_0058);
    const ccomplex_t IT_0226 = (-6)*conj(IT_0222);
    const ccomplex_t IT_0227 = IT_0225 + IT_0226;
    const ccomplex_t IT_0228 = IT_0072*IT_0149;
    const ccomplex_t IT_0229 = IT_0082*IT_0207*IT_0228;
    const ccomplex_t IT_0230 = IT_0092*IT_0163;
    const ccomplex_t IT_0231 = IT_0102*IT_0210*IT_0230;
    const ccomplex_t IT_0232 = IT_0132*IT_0156;
    const ccomplex_t IT_0233 = IT_0142*IT_0216*IT_0232;
    const ccomplex_t IT_0234 = m_s*IT_0161;
    const ccomplex_t IT_0235 = IT_0121*IT_0172;
    const ccomplex_t IT_0236 = IT_0105*IT_0213*IT_0235;
    const ccomplex_t IT_0237 = m_s*IT_0154;
    const ccomplex_t IT_0238 = m_s*IT_0168;
    const ccomplex_t IT_0239 = m_s*IT_0175;
    const ccomplex_t IT_0240 = -IT_0229 + -IT_0231 + -IT_0233 + IT_0234 + 
      -IT_0236 + IT_0237 + IT_0238 + IT_0239;
    const ccomplex_t IT_0241 = (-18)*IT_0222;
    const ccomplex_t IT_0242 = (-6)*IT_0058;
    const ccomplex_t IT_0243 = IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = (-18)*conj(IT_0222);
    const ccomplex_t IT_0245 = (-6)*conj(IT_0058);
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = 6*s_34;
    const ccomplex_t IT_0248 = -IT_0008;
    const ccomplex_t IT_0249 = s_34 + IT_0248;
    const ccomplex_t IT_0250 = 6*IT_0027;
    const ccomplex_t IT_0251 = 6*conj(IT_0027);
    const ccomplex_t IT_0252 = (-6)*IT_0194;
    const ccomplex_t IT_0253 = IT_0178 + IT_0180 + IT_0252;
    const ccomplex_t IT_0254 = conj(IT_0064) + conj(IT_0147);
    const ccomplex_t IT_0255 = conj(IT_0145) + conj(IT_0176);
    const ccomplex_t IT_0256 = IT_0145 + IT_0176;
    const ccomplex_t IT_0257 = -IT_0204;
    const ccomplex_t IT_0258 = IT_0064 + IT_0147;
    const ccomplex_t IT_0259 = -IT_0191;
    const ccomplex_t IT_0260 = (-18)*IT_0194;
    const ccomplex_t IT_0261 = IT_0182 + IT_0183 + IT_0260;
    const ccomplex_t IT_0262 = m_s*s_13;
    const ccomplex_t IT_0263 = 6*conj(IT_0222);
    const ccomplex_t IT_0264 = 18*conj(IT_0058);
    const ccomplex_t IT_0265 = 6*conj(IT_0058);
    const ccomplex_t IT_0266 = 18*conj(IT_0222);
    const ccomplex_t IT_0267 = 18*conj(IT_0240);
    const ccomplex_t IT_0268 = 6*conj(IT_0240);
    const ccomplex_t IT_0269 = 18*conj(IT_0059);
    const ccomplex_t IT_0270 = 0.166666666666667*IT_0058*(18*conj(IT_0064) + 6
      *conj(IT_0145)) + 0.166666666666667*IT_0059*(18*conj(IT_0147) + 6*conj
      (IT_0176)) + (conj(IT_0064) + 3*conj(IT_0145))*IT_0222 + (conj(IT_0147) +
       3*conj(IT_0176))*IT_0240 + 0.166666666666667*IT_0064*(IT_0263 + IT_0264) 
      + 0.166666666666667*IT_0145*(IT_0265 + IT_0266) + IT_0176*(conj(IT_0059) +
       0.166666666666667*IT_0267) + 0.166666666666667*IT_0147*(IT_0268 + IT_0269);
    const ccomplex_t IT_0271 = -s_14;
    const ccomplex_t IT_0272 = s_13 + IT_0271;
    const ccomplex_t IT_0273 = m_s*IT_0272;
    const ccomplex_t IT_0274 = 18*IT_0027;
    const ccomplex_t IT_0275 = 18*conj(IT_0027);
    const ccomplex_t IT_0276 = (-18)*IT_0240;
    const ccomplex_t IT_0277 = (-6)*IT_0240;
    const ccomplex_t IT_0278 = IT_0250*IT_0255 + IT_0251*IT_0256 + IT_0254
      *IT_0274 + IT_0258*IT_0275 + IT_0191*((-12)*IT_0027 + 12*conj(IT_0027) + (
      -6)*IT_0059 + 6*conj(IT_0059) + IT_0241 + IT_0242 + IT_0265 + IT_0266 +
       IT_0267 + IT_0276) + IT_0204*((-36)*IT_0027 + 36*conj(IT_0027) + (-18)
      *IT_0059 + IT_0205 + IT_0223 + IT_0263 + IT_0264 + IT_0268 + IT_0269 +
       IT_0277);
    const ccomplex_t IT_0279 = m_s*s_14;
    const ccomplex_t IT_0280 = (-18)*IT_0064;
    const ccomplex_t IT_0281 = (-6)*IT_0145;
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = (-18)*conj(IT_0064);
    const ccomplex_t IT_0284 = (-6)*conj(IT_0145);
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = (-0.166666666666667)*IT_0147;
    const ccomplex_t IT_0287 = (-0.166666666666667)*conj(IT_0147);
    const ccomplex_t IT_0288 = (IT_0064 + 3*IT_0145)*conj(IT_0240) + (
      -0.166666666666667)*conj(IT_0176)*IT_0243 + (-0.166666666666667)*IT_0176
      *IT_0246 + (-0.166666666666667)*conj(IT_0145)*IT_0276 + (
      -0.166666666666667)*conj(IT_0064)*IT_0277 + (-0.166666666666667)*conj
      (IT_0059)*IT_0282 + (-0.166666666666667)*IT_0059*IT_0285 + IT_0227*IT_0286
       + IT_0224*IT_0287;
    const ccomplex_t IT_0289 = (36*s_34 + (-36)*IT_0008)*IT_0027*conj(IT_0027)
       + (conj(IT_0027)*(IT_0058 + IT_0059) + IT_0027*(conj(IT_0058) + conj
      (IT_0059)))*IT_0062 + (conj(IT_0064)*IT_0145 + IT_0064*conj(IT_0145) +
       conj(IT_0147)*IT_0176 + IT_0147*conj(IT_0176))*IT_0181 + (IT_0064*conj
      (IT_0064) + IT_0145*conj(IT_0145) + IT_0147*conj(IT_0147) + IT_0176*conj
      (IT_0176))*IT_0184 + -IT_0192*IT_0197 + (-2)*IT_0204*(IT_0191*((-12)
      *IT_0177 + 24*IT_0179 + (-12)*IT_0194) + 0.5*IT_0197*IT_0204) + IT_0061*
      (IT_0058*conj(IT_0058) + IT_0059*conj(IT_0059) + IT_0222*conj(IT_0222) +
       IT_0240*conj(IT_0240)) + IT_0008*(conj(IT_0059)*IT_0224 + IT_0059*IT_0227
       + conj(IT_0240)*IT_0243 + IT_0240*IT_0246) + (conj(IT_0058)*IT_0222 +
       IT_0058*conj(IT_0222) + conj(IT_0059)*IT_0240 + IT_0059*conj(IT_0240))
      *IT_0247 + IT_0249*((conj(IT_0222) + conj(IT_0240))*IT_0250 + (IT_0222 +
       IT_0240)*IT_0251) + IT_0253*(IT_0191*IT_0254 + IT_0204*IT_0255 + IT_0256
      *IT_0257 + IT_0258*IT_0259) + (IT_0204*IT_0254 + IT_0191*IT_0255 + IT_0257
      *IT_0258 + IT_0256*IT_0259)*IT_0261 + 6*IT_0262*IT_0270 + IT_0273*IT_0278 
      + (-6)*IT_0194*((conj(IT_0064) + 3*conj(IT_0145))*IT_0176 + (IT_0064 + 3
      *IT_0145)*conj(IT_0176) + IT_0285*IT_0286 + IT_0282*IT_0287) + (-6)
      *IT_0279*IT_0288;
    return create_ccomplex_return(IT_0289);
}

