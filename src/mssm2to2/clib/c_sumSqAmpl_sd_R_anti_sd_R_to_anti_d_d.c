#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_anti_sd_R_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_anti_sd_R_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
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
    const creal_t m_sd_R = param->m_sd_R;
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
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = pow(m_d, 2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0001*IT_0003 + -IT_0010*IT_0011) + (-3)*IT_0000*IT_0002*IT_0003
      *IT_0005*IT_0012);
    const ccomplex_t IT_0014 = 0.333333333333333*IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_sd_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0011;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0003*IT_0010 + IT_0001*IT_0011) + (-3)*IT_0000*IT_0002*IT_0005*IT_0011
      *IT_0012);
    const ccomplex_t IT_0023 = (-0.333333333333333)*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0019 + -IT_0027;
    const ccomplex_t IT_0029 = cpow(IT_0008, -1);
    const ccomplex_t IT_0030 = IT_0004*IT_0029;
    const ccomplex_t IT_0031 = e_em*IT_0030;
    const ccomplex_t IT_0032 = IT_0005*IT_0008;
    const ccomplex_t IT_0033 = e_em*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0031 + 3*IT_0033);
    const ccomplex_t IT_0035 = (-0.166666666666667)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0029;
    const ccomplex_t IT_0040 = (-0.666666666666667)*IT_0039;
    const ccomplex_t IT_0041 = m_d*IT_0040;
    const ccomplex_t IT_0042 = (-0.333333333333333)*IT_0039;
    const ccomplex_t IT_0043 = m_d*IT_0042;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = IT_0041 + IT_0044;
    const ccomplex_t IT_0046 = IT_0038*IT_0045;
    const ccomplex_t IT_0047 = 0.333333333333333*IT_0039;
    const ccomplex_t IT_0048 = IT_0036*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0045*IT_0049;
    const ccomplex_t IT_0051 = pow(m_sG, 2);
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0012 + IT_0016 + -IT_0051 +
       reg_prop, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0054 = cpow(IT_0053, 2);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0052*IT_0055;
    const ccomplex_t IT_0057 = m_d*IT_0056;
    const ccomplex_t IT_0058 = -IT_0046 + IT_0050 + 0.5*IT_0057;
    const ccomplex_t IT_0059 = (IT_0041 + -IT_0043)*(IT_0038 + -IT_0049);
    const ccomplex_t IT_0060 = (-18)*IT_0012;
    const ccomplex_t IT_0061 = 18*s_34;
    const ccomplex_t IT_0062 = IT_0060 + IT_0061;
    const ccomplex_t IT_0063 = IT_0038*IT_0040;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0012 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0012 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0012 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_23 + IT_0012 + IT_0016 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = -IT_0072 + -IT_0080 + -IT_0088 + -IT_0096;
    const ccomplex_t IT_0098 = IT_0040*IT_0049;
    const ccomplex_t IT_0099 = (-0.5)*IT_0056 + -IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0029;
    const ccomplex_t IT_0101 = (-0.333333333333333)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0029;
    const ccomplex_t IT_0103 = (-0.333333333333333)*IT_0102;
    const ccomplex_t IT_0104 = IT_0101*IT_0103;
    const ccomplex_t IT_0105 = IT_0086*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0029;
    const ccomplex_t IT_0108 = (-0.333333333333333)*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0029;
    const ccomplex_t IT_0110 = (-0.333333333333333)*IT_0109;
    const ccomplex_t IT_0111 = IT_0108*IT_0110;
    const ccomplex_t IT_0112 = IT_0094*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0029;
    const ccomplex_t IT_0115 = (-0.333333333333333)*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0029;
    const ccomplex_t IT_0117 = (-0.333333333333333)*IT_0116;
    const ccomplex_t IT_0118 = IT_0115*IT_0117;
    const ccomplex_t IT_0119 = IT_0070*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0029;
    const ccomplex_t IT_0122 = (-0.333333333333333)*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0029;
    const ccomplex_t IT_0124 = (-0.333333333333333)*IT_0123;
    const ccomplex_t IT_0125 = IT_0122*IT_0124;
    const ccomplex_t IT_0126 = IT_0078*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = 0.166666666666667*IT_0056 + -IT_0106 + -IT_0113
       + -IT_0120 + -IT_0127;
    const ccomplex_t IT_0129 = s_34*IT_0016;
    const ccomplex_t IT_0130 = (-6)*IT_0129;
    const ccomplex_t IT_0131 = s_13*s_14;
    const ccomplex_t IT_0132 = 12*IT_0131;
    const ccomplex_t IT_0133 = IT_0130 + IT_0132;
    const ccomplex_t IT_0134 = (-18)*IT_0129;
    const ccomplex_t IT_0135 = 36*IT_0131;
    const ccomplex_t IT_0136 = IT_0134 + IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0138 = (-2)*IT_0137;
    const ccomplex_t IT_0139 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0140 = IT_0137*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = IT_0138*IT_0141;
    const ccomplex_t IT_0143 = 0.25*IT_0142;
    const ccomplex_t IT_0144 = cpow(IT_0143, 2);
    const ccomplex_t IT_0145 = (-36)*IT_0129;
    const ccomplex_t IT_0146 = IT_0012*IT_0016;
    const ccomplex_t IT_0147 = (-36)*IT_0146;
    const ccomplex_t IT_0148 = 72*IT_0131;
    const ccomplex_t IT_0149 = IT_0145 + IT_0147 + IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0151 = 0.666666666666667*IT_0150;
    const ccomplex_t IT_0152 = (-0.333333333333333)*IT_0150;
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = IT_0139*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = (-0.0833333333333333)*IT_0142 + 0.5*IT_0155;
    const ccomplex_t IT_0157 = 3*conj(IT_0128);
    const ccomplex_t IT_0158 = conj(IT_0099) + IT_0157;
    const ccomplex_t IT_0159 = IT_0068*IT_0117;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0161 = IT_0070*IT_0159*IT_0160;
    const ccomplex_t IT_0162 = IT_0076*IT_0124;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0164 = IT_0078*IT_0162*IT_0163;
    const ccomplex_t IT_0165 = IT_0084*IT_0103;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0167 = IT_0086*IT_0165*IT_0166;
    const ccomplex_t IT_0168 = IT_0092*IT_0110;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0170 = IT_0094*IT_0168*IT_0169;
    const ccomplex_t IT_0171 = m_d*IT_0072;
    const ccomplex_t IT_0172 = m_d*IT_0080;
    const ccomplex_t IT_0173 = m_d*IT_0088;
    const ccomplex_t IT_0174 = m_d*IT_0096;
    const ccomplex_t IT_0175 = -IT_0161 + -IT_0164 + -IT_0167 + -IT_0170 +
       IT_0171 + IT_0172 + IT_0173 + IT_0174;
    const ccomplex_t IT_0176 = IT_0066*IT_0115;
    const ccomplex_t IT_0177 = IT_0070*IT_0160*IT_0176;
    const ccomplex_t IT_0178 = IT_0074*IT_0122;
    const ccomplex_t IT_0179 = IT_0078*IT_0163*IT_0178;
    const ccomplex_t IT_0180 = IT_0082*IT_0101;
    const ccomplex_t IT_0181 = IT_0086*IT_0166*IT_0180;
    const ccomplex_t IT_0182 = IT_0090*IT_0108;
    const ccomplex_t IT_0183 = IT_0094*IT_0169*IT_0182;
    const ccomplex_t IT_0184 = m_d*IT_0106;
    const ccomplex_t IT_0185 = m_d*IT_0113;
    const ccomplex_t IT_0186 = m_d*IT_0120;
    const ccomplex_t IT_0187 = m_d*IT_0127;
    const ccomplex_t IT_0188 = (-0.166666666666667)*IT_0057 + -IT_0177 + 
      -IT_0179 + -IT_0181 + -IT_0183 + IT_0184 + IT_0185 + IT_0186 + IT_0187;
    const ccomplex_t IT_0189 = 6*s_34;
    const ccomplex_t IT_0190 = -IT_0012;
    const ccomplex_t IT_0191 = s_34 + IT_0190;
    const ccomplex_t IT_0192 = 6*IT_0028;
    const ccomplex_t IT_0193 = 6*conj(IT_0028);
    const ccomplex_t IT_0194 = m_d*s_13;
    const ccomplex_t IT_0195 = 6*conj(IT_0175);
    const ccomplex_t IT_0196 = 18*conj(IT_0175);
    const ccomplex_t IT_0197 = 6*conj(IT_0188);
    const ccomplex_t IT_0198 = 18*conj(IT_0058);
    const ccomplex_t IT_0199 = 6*conj(IT_0058);
    const ccomplex_t IT_0200 = 18*conj(IT_0188);
    const ccomplex_t IT_0201 = 3*IT_0059*(conj(IT_0064) + 0.333333333333333
      *conj(IT_0097)) + 3*IT_0058*(conj(IT_0099) + 0.333333333333333*conj
      (IT_0128)) + (conj(IT_0064) + 3*conj(IT_0097))*IT_0175 + IT_0158*IT_0188 +
       3*IT_0064*(conj(IT_0059) + 0.0555555555555556*IT_0195) + IT_0097*(conj
      (IT_0059) + 0.166666666666667*IT_0196) + 0.166666666666667*IT_0099*
      (IT_0197 + IT_0198) + 0.166666666666667*IT_0128*(IT_0199 + IT_0200);
    const ccomplex_t IT_0202 = (-18)*IT_0058;
    const ccomplex_t IT_0203 = (-6)*IT_0188;
    const ccomplex_t IT_0204 = IT_0202 + IT_0203;
    const ccomplex_t IT_0205 = (-18)*conj(IT_0058);
    const ccomplex_t IT_0206 = (-6)*conj(IT_0188);
    const ccomplex_t IT_0207 = IT_0205 + IT_0206;
    const ccomplex_t IT_0208 = (-18)*IT_0175;
    const ccomplex_t IT_0209 = (-18)*conj(IT_0175);
    const ccomplex_t IT_0210 = (-6)*IT_0175;
    const ccomplex_t IT_0211 = (-6)*conj(IT_0175);
    const ccomplex_t IT_0212 = m_d*s_14;
    const ccomplex_t IT_0213 = (-18)*IT_0188;
    const ccomplex_t IT_0214 = (-6)*IT_0058;
    const ccomplex_t IT_0215 = IT_0097*(conj(IT_0058) + 3*conj(IT_0188)) + (
      -0.166666666666667)*conj(IT_0064)*IT_0204 + (-0.166666666666667)*IT_0064
      *IT_0207 + conj(IT_0128)*(IT_0059 + (-0.166666666666667)*IT_0208) +
       IT_0128*(conj(IT_0059) + (-0.166666666666667)*IT_0209) + 3*conj(IT_0099)*
      (IT_0059 + (-0.0555555555555556)*IT_0210) + 3*IT_0099*(conj(IT_0059) + (
      -0.0555555555555556)*IT_0211) + (-0.166666666666667)*conj(IT_0097)*
      (IT_0213 + IT_0214);
    const ccomplex_t IT_0216 = -s_14;
    const ccomplex_t IT_0217 = s_13 + IT_0216;
    const ccomplex_t IT_0218 = m_d*IT_0217;
    const ccomplex_t IT_0219 = conj(IT_0097) + conj(IT_0128);
    const ccomplex_t IT_0220 = IT_0097 + IT_0128;
    const ccomplex_t IT_0221 = conj(IT_0064) + conj(IT_0099);
    const ccomplex_t IT_0222 = 18*IT_0028;
    const ccomplex_t IT_0223 = IT_0064 + IT_0099;
    const ccomplex_t IT_0224 = 18*conj(IT_0028);
    const ccomplex_t IT_0225 = IT_0156*((-36)*IT_0028 + 36*conj(IT_0028) + (
      -18)*IT_0059 + 18*conj(IT_0059) + IT_0195 + IT_0197 + IT_0198 + IT_0202 +
       IT_0203 + IT_0210) + IT_0143*((-12)*IT_0028 + 12*conj(IT_0028) + (-6)
      *IT_0059 + 6*conj(IT_0059) + IT_0196 + IT_0199 + IT_0200 + IT_0208 +
       IT_0213 + IT_0214) + IT_0192*IT_0219 + IT_0193*IT_0220 + IT_0221*IT_0222 
      + IT_0223*IT_0224;
    const ccomplex_t IT_0226 = (-6)*IT_0146;
    const ccomplex_t IT_0227 = IT_0130 + IT_0132 + IT_0226;
    const ccomplex_t IT_0228 = -IT_0156;
    const ccomplex_t IT_0229 = -IT_0143;
    const ccomplex_t IT_0230 = (-18)*IT_0146;
    const ccomplex_t IT_0231 = IT_0134 + IT_0135 + IT_0230;
    const ccomplex_t IT_0232 = (36*s_34 + (-36)*IT_0012)*IT_0028*conj(IT_0028)
       + (conj(IT_0028)*(IT_0058 + IT_0059) + IT_0028*(conj(IT_0058) + conj
      (IT_0059)))*IT_0062 + (conj(IT_0064)*IT_0097 + IT_0064*conj(IT_0097) +
       conj(IT_0099)*IT_0128 + IT_0099*conj(IT_0128))*IT_0133 + (IT_0064*conj
      (IT_0064) + IT_0097*conj(IT_0097) + IT_0099*conj(IT_0099) + IT_0128*conj
      (IT_0128))*IT_0136 + -IT_0144*IT_0149 + (-2)*IT_0156*(IT_0143*((-12)
      *IT_0129 + 24*IT_0131 + (-12)*IT_0146) + 0.5*IT_0149*IT_0156) + (-6)
      *IT_0146*(conj(IT_0097)*(IT_0099 + 3*IT_0128) + (-0.166666666666667)*conj
      (IT_0064)*((-18)*IT_0099 + (-6)*IT_0128) + (-0.166666666666667)*IT_0064*((
      -18)*conj(IT_0099) + (-6)*conj(IT_0128)) + IT_0097*IT_0158) + IT_0061*
      (IT_0058*conj(IT_0058) + IT_0059*conj(IT_0059) + IT_0175*conj(IT_0175) +
       IT_0188*conj(IT_0188)) + (conj(IT_0059)*IT_0175 + IT_0059*conj(IT_0175) +
       conj(IT_0058)*IT_0188 + IT_0058*conj(IT_0188))*IT_0189 + IT_0191*((conj
      (IT_0175) + conj(IT_0188))*IT_0192 + (IT_0175 + IT_0188)*IT_0193) + 6
      *IT_0194*IT_0201 + IT_0012*(conj(IT_0059)*IT_0204 + IT_0059*IT_0207 + conj
      (IT_0188)*IT_0208 + IT_0188*IT_0209 + conj(IT_0058)*IT_0210 + IT_0058
      *IT_0211) + (-6)*IT_0212*IT_0215 + IT_0218*IT_0225 + IT_0227*(IT_0156
      *IT_0219 + IT_0143*IT_0221 + IT_0220*IT_0228 + IT_0223*IT_0229) + (IT_0143
      *IT_0219 + IT_0156*IT_0221 + IT_0223*IT_0228 + IT_0220*IT_0229)*IT_0231;
    return create_ccomplex_return(IT_0232);
}

