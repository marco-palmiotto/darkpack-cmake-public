#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_anti_ss_R_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_anti_ss_R_to_anti_s_s(
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
    const creal_t m_ss_R = param->m_ss_R;
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
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = cos(alpha);
    const ccomplex_t IT_0011 = sin(beta);
    const ccomplex_t IT_0012 = pow(m_s, 2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0001*IT_0003 + IT_0010*IT_0011) + (-3)*IT_0000*IT_0002*IT_0003*IT_0005
      *IT_0012);
    const ccomplex_t IT_0014 = (-0.333333333333333)*IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_ss_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0005*IT_0010;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0001*IT_0010 + -IT_0003*IT_0011) + (-3)*IT_0000*IT_0002*IT_0005
      *IT_0010*IT_0012);
    const ccomplex_t IT_0023 = 0.333333333333333*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0019 + -IT_0027;
    const ccomplex_t IT_0029 = pow(m_sG, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0012 + IT_0016 + -IT_0029 +
       reg_prop, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0032 = cpow(IT_0031, 2);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = m_s*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0008, -1);
    const ccomplex_t IT_0037 = IT_0004*IT_0036;
    const ccomplex_t IT_0038 = e_em*IT_0037;
    const ccomplex_t IT_0039 = IT_0005*IT_0008;
    const ccomplex_t IT_0040 = e_em*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0038 + 3*IT_0040);
    const ccomplex_t IT_0042 = (-0.166666666666667)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0036;
    const ccomplex_t IT_0047 = (-0.666666666666667)*IT_0046;
    const ccomplex_t IT_0048 = m_s*IT_0047;
    const ccomplex_t IT_0049 = (-0.333333333333333)*IT_0046;
    const ccomplex_t IT_0050 = m_s*IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0048 + IT_0051;
    const ccomplex_t IT_0053 = IT_0045*IT_0052;
    const ccomplex_t IT_0054 = 0.333333333333333*IT_0046;
    const ccomplex_t IT_0055 = IT_0043*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = IT_0052*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0035 + -IT_0053 + IT_0057;
    const ccomplex_t IT_0059 = (IT_0048 + -IT_0050)*(IT_0045 + -IT_0056);
    const ccomplex_t IT_0060 = (-18)*IT_0012;
    const ccomplex_t IT_0061 = 18*s_34;
    const ccomplex_t IT_0062 = IT_0060 + IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0064 = (-2)*IT_0063;
    const ccomplex_t IT_0065 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0064*IT_0067;
    const ccomplex_t IT_0069 = 0.25*IT_0068;
    const ccomplex_t IT_0070 = cpow(IT_0069, 2);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0072 = 0.666666666666667*IT_0071;
    const ccomplex_t IT_0073 = (-0.333333333333333)*IT_0071;
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = IT_0065*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = (-0.0833333333333333)*IT_0068 + 0.5*IT_0076;
    const ccomplex_t IT_0078 = cpow(IT_0077, 2);
    const ccomplex_t IT_0079 = s_34*IT_0016;
    const ccomplex_t IT_0080 = (-36)*IT_0079;
    const ccomplex_t IT_0081 = IT_0012*IT_0016;
    const ccomplex_t IT_0082 = (-36)*IT_0081;
    const ccomplex_t IT_0083 = s_13*s_14;
    const ccomplex_t IT_0084 = 72*IT_0083;
    const ccomplex_t IT_0085 = IT_0080 + IT_0082 + IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0036;
    const ccomplex_t IT_0087 = (-0.333333333333333)*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0092 = cpow((-2)*s_23 + IT_0012 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0090*IT_0091*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0036;
    const ccomplex_t IT_0095 = (-0.333333333333333)*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = IT_0095*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0100 = cpow((-2)*s_23 + IT_0012 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0101 = IT_0098*IT_0099*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0036;
    const ccomplex_t IT_0103 = (-0.333333333333333)*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = IT_0103*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0108 = cpow((-2)*s_23 + IT_0012 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0106*IT_0107*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0036;
    const ccomplex_t IT_0111 = (-0.333333333333333)*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0116 = cpow((-2)*s_23 + IT_0012 + IT_0016 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0114*IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0097*IT_0119;
    const ccomplex_t IT_0121 = IT_0100*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = m_s*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = IT_0105*IT_0125;
    const ccomplex_t IT_0127 = IT_0108*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = m_s*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = IT_0113*IT_0131;
    const ccomplex_t IT_0133 = IT_0116*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = m_s*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0089*IT_0137;
    const ccomplex_t IT_0139 = IT_0092*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = m_s*IT_0140;
    const ccomplex_t IT_0142 = -IT_0093 + -IT_0101 + -IT_0109 + -IT_0117 +
       IT_0123 + IT_0129 + IT_0135 + IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0036;
    const ccomplex_t IT_0144 = (-0.333333333333333)*IT_0143;
    const ccomplex_t IT_0145 = IT_0119*IT_0144;
    const ccomplex_t IT_0146 = IT_0099*IT_0100*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0036;
    const ccomplex_t IT_0148 = (-0.333333333333333)*IT_0147;
    const ccomplex_t IT_0149 = IT_0125*IT_0148;
    const ccomplex_t IT_0150 = IT_0107*IT_0108*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0036;
    const ccomplex_t IT_0152 = (-0.333333333333333)*IT_0151;
    const ccomplex_t IT_0153 = IT_0131*IT_0152;
    const ccomplex_t IT_0154 = IT_0115*IT_0116*IT_0153;
    const ccomplex_t IT_0155 = IT_0103*IT_0148;
    const ccomplex_t IT_0156 = IT_0108*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = m_s*IT_0157;
    const ccomplex_t IT_0159 = IT_0111*IT_0152;
    const ccomplex_t IT_0160 = IT_0116*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = m_s*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0036;
    const ccomplex_t IT_0164 = (-0.333333333333333)*IT_0163;
    const ccomplex_t IT_0165 = IT_0137*IT_0164;
    const ccomplex_t IT_0166 = IT_0091*IT_0092*IT_0165;
    const ccomplex_t IT_0167 = IT_0087*IT_0164;
    const ccomplex_t IT_0168 = IT_0092*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = m_s*IT_0169;
    const ccomplex_t IT_0171 = IT_0095*IT_0144;
    const ccomplex_t IT_0172 = IT_0100*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = m_s*IT_0173;
    const ccomplex_t IT_0175 = (-0.166666666666667)*IT_0035 + -IT_0146 + 
      -IT_0150 + -IT_0154 + IT_0158 + IT_0162 + -IT_0166 + IT_0170 + IT_0174;
    const ccomplex_t IT_0176 = 6*s_34;
    const ccomplex_t IT_0177 = -IT_0012;
    const ccomplex_t IT_0178 = s_34 + IT_0177;
    const ccomplex_t IT_0179 = 6*IT_0028;
    const ccomplex_t IT_0180 = 6*conj(IT_0028);
    const ccomplex_t IT_0181 = (-18)*IT_0058;
    const ccomplex_t IT_0182 = (-6)*IT_0175;
    const ccomplex_t IT_0183 = IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = (-18)*conj(IT_0058);
    const ccomplex_t IT_0185 = (-6)*conj(IT_0175);
    const ccomplex_t IT_0186 = IT_0184 + IT_0185;
    const ccomplex_t IT_0187 = (-18)*IT_0142;
    const ccomplex_t IT_0188 = (-18)*conj(IT_0142);
    const ccomplex_t IT_0189 = (-6)*IT_0142;
    const ccomplex_t IT_0190 = (-6)*conj(IT_0142);
    const ccomplex_t IT_0191 = (-2)*IT_0069;
    const ccomplex_t IT_0192 = (-18)*IT_0079;
    const ccomplex_t IT_0193 = (-18)*IT_0081;
    const ccomplex_t IT_0194 = 36*IT_0083;
    const ccomplex_t IT_0195 = IT_0192 + IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = IT_0045*IT_0047;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = IT_0047*IT_0056;
    const ccomplex_t IT_0199 = -IT_0198;
    const ccomplex_t IT_0200 = conj(IT_0197) + conj(IT_0199);
    const ccomplex_t IT_0201 = -IT_0122 + -IT_0128 + -IT_0134 + -IT_0140;
    const ccomplex_t IT_0202 = -IT_0157 + -IT_0161 + -IT_0169 + -IT_0173;
    const ccomplex_t IT_0203 = conj(IT_0201) + conj(IT_0202);
    const ccomplex_t IT_0204 = 0.166666666666667*IT_0034;
    const ccomplex_t IT_0205 = (-0.5)*IT_0034;
    const ccomplex_t IT_0206 = (-2)*IT_0077;
    const ccomplex_t IT_0207 = IT_0201 + IT_0202;
    const ccomplex_t IT_0208 = -IT_0069;
    const ccomplex_t IT_0209 = IT_0197 + IT_0199;
    const ccomplex_t IT_0210 = -IT_0077;
    const ccomplex_t IT_0211 = (-6)*IT_0079;
    const ccomplex_t IT_0212 = (-6)*IT_0081;
    const ccomplex_t IT_0213 = 12*IT_0083;
    const ccomplex_t IT_0214 = IT_0211 + IT_0212 + IT_0213;
    const ccomplex_t IT_0215 = -s_14;
    const ccomplex_t IT_0216 = s_13 + IT_0215;
    const ccomplex_t IT_0217 = m_s*IT_0216;
    const ccomplex_t IT_0218 = 18*IT_0028;
    const ccomplex_t IT_0219 = 18*conj(IT_0028);
    const ccomplex_t IT_0220 = 6*conj(IT_0058);
    const ccomplex_t IT_0221 = 18*conj(IT_0142);
    const ccomplex_t IT_0222 = 18*conj(IT_0175);
    const ccomplex_t IT_0223 = (-18)*IT_0175;
    const ccomplex_t IT_0224 = (-6)*IT_0058;
    const ccomplex_t IT_0225 = 6*conj(IT_0142);
    const ccomplex_t IT_0226 = 6*conj(IT_0175);
    const ccomplex_t IT_0227 = 18*conj(IT_0058);
    const ccomplex_t IT_0228 = (-18)*IT_0059;
    const ccomplex_t IT_0229 = (-0.166666666666667)*IT_0179*IT_0203 + (
      -0.166666666666667)*IT_0180*(IT_0201 + IT_0202 + IT_0204) + IT_0028*
      (IT_0204 + 3*IT_0205) + (-0.166666666666667)*IT_0200*IT_0218 + (
      -0.166666666666667)*(IT_0205 + IT_0209)*IT_0219 + (-0.166666666666667)
      *IT_0069*((-12)*IT_0028 + 12*conj(IT_0028) + (-6)*IT_0059 + 6*conj(IT_0059
      ) + IT_0187 + IT_0220 + IT_0221 + IT_0222 + IT_0223 + IT_0224) + (
      -0.166666666666667)*IT_0077*((-36)*IT_0028 + 36*conj(IT_0028) + 18*conj
      (IT_0059) + IT_0181 + IT_0182 + IT_0189 + IT_0225 + IT_0226 + IT_0227 +
       IT_0228);
    const ccomplex_t IT_0230 = IT_0192 + IT_0194;
    const ccomplex_t IT_0231 = -IT_0202;
    const ccomplex_t IT_0232 = -IT_0204;
    const ccomplex_t IT_0233 = -conj(IT_0199);
    const ccomplex_t IT_0234 = IT_0211 + IT_0213;
    const ccomplex_t IT_0235 = m_s*s_13;
    const ccomplex_t IT_0236 = conj(IT_0199) + 3*conj(IT_0202);
    const ccomplex_t IT_0237 = IT_0220 + IT_0222;
    const ccomplex_t IT_0238 = IT_0226 + IT_0227;
    const ccomplex_t IT_0239 = 6*IT_0201;
    const ccomplex_t IT_0240 = 18*IT_0197;
    const ccomplex_t IT_0241 = IT_0239 + IT_0240;
    const ccomplex_t IT_0242 = 6*IT_0142;
    const ccomplex_t IT_0243 = 18*IT_0059;
    const ccomplex_t IT_0244 = 18*IT_0142;
    const ccomplex_t IT_0245 = 0.166666666666667*IT_0244;
    const ccomplex_t IT_0246 = IT_0058*(18*conj(IT_0199) + 6*conj(IT_0202)) +
       IT_0201*IT_0221 + IT_0197*IT_0225 + 6*IT_0175*IT_0236 + IT_0204*(IT_0223 
      + IT_0224 + IT_0237) + IT_0202*IT_0237 + IT_0205*(IT_0181 + IT_0182 +
       IT_0238) + IT_0199*IT_0238 + conj(IT_0059)*IT_0241 + conj(IT_0197)*
      (IT_0242 + IT_0243) + 6*conj(IT_0201)*(IT_0059 + IT_0245);
    const ccomplex_t IT_0247 = (-18)*conj(IT_0201);
    const ccomplex_t IT_0248 = (-6)*conj(IT_0197);
    const ccomplex_t IT_0249 = (-0.166666666666667)*IT_0197;
    const ccomplex_t IT_0250 = (-18)*conj(IT_0197);
    const ccomplex_t IT_0251 = (-6)*conj(IT_0201);
    const ccomplex_t IT_0252 = (-0.166666666666667)*IT_0199;
    const ccomplex_t IT_0253 = (-0.166666666666667)*IT_0205;
    const ccomplex_t IT_0254 = IT_0201*IT_0236 + (-0.166666666666667)*IT_0202*
      (IT_0247 + IT_0248) + IT_0232*(IT_0197 + 3*IT_0201 + 0.166666666666667
      *IT_0247 + 0.166666666666667*IT_0248) + ((-18)*conj(IT_0199) + (-6)*conj
      (IT_0202))*IT_0249 + (IT_0250 + IT_0251)*IT_0252 + (IT_0241 + IT_0250 +
       IT_0251)*IT_0253;
    const ccomplex_t IT_0255 = m_s*s_14;
    const ccomplex_t IT_0256 = (-18)*conj(IT_0059);
    const ccomplex_t IT_0257 = (-0.166666666666667)*IT_0183*conj(IT_0197) + 
      (conj(IT_0058) + 3*conj(IT_0175))*IT_0201 + (conj(IT_0059) + (
      -0.166666666666667)*IT_0188)*IT_0202 + (IT_0059 + (-0.166666666666667)
      *IT_0187)*conj(IT_0202) + (-0.166666666666667)*conj(IT_0201)*(IT_0223 +
       IT_0224) + (-0.166666666666667)*conj(IT_0199)*(IT_0189 + IT_0228) +
       IT_0232*(IT_0059 + -conj(IT_0059) + 0.166666666666667*IT_0188 + IT_0245) 
      + IT_0186*IT_0249 + IT_0252*(IT_0190 + IT_0256) + IT_0253*(IT_0190 +
       IT_0242 + IT_0243 + IT_0256);
    const ccomplex_t IT_0258 = (36*s_34 + (-36)*IT_0012)*IT_0028*conj(IT_0028)
       + (conj(IT_0028)*(IT_0058 + IT_0059) + IT_0028*(conj(IT_0058) + conj
      (IT_0059)))*IT_0062 + -(IT_0070 + IT_0078)*IT_0085 + IT_0061*(IT_0058*conj
      (IT_0058) + IT_0059*conj(IT_0059) + IT_0142*conj(IT_0142) + IT_0175*conj
      (IT_0175)) + (conj(IT_0059)*IT_0142 + IT_0059*conj(IT_0142) + conj(IT_0058
      )*IT_0175 + IT_0058*conj(IT_0175))*IT_0176 + IT_0178*((conj(IT_0142) +
       conj(IT_0175))*IT_0179 + (IT_0142 + IT_0175)*IT_0180) + IT_0012*(conj
      (IT_0059)*IT_0183 + IT_0059*IT_0186 + conj(IT_0175)*IT_0187 + IT_0175
      *IT_0188 + conj(IT_0058)*IT_0189 + IT_0058*IT_0190) + IT_0077*((-12)
      *IT_0079 + (-12)*IT_0081 + 24*IT_0083)*IT_0191 + IT_0195*(IT_0077*IT_0200 
      + IT_0069*IT_0203 + IT_0191*IT_0204 + IT_0205*IT_0206 + IT_0207*IT_0208 +
       IT_0209*IT_0210) + (IT_0069*IT_0200 + IT_0077*IT_0203 + IT_0191*IT_0205 +
       IT_0204*IT_0206 + IT_0208*IT_0209 + IT_0207*IT_0210)*IT_0214 + (-6)
      *IT_0217*IT_0229 + IT_0230*(IT_0197*conj(IT_0197) + IT_0199*conj(IT_0199) 
      + IT_0201*conj(IT_0201) + IT_0202*conj(IT_0202) + -(conj(IT_0202) + 
      -IT_0204 + IT_0231)*IT_0232 + -IT_0205*(IT_0199 + IT_0205 + IT_0233)) + 
      (conj(IT_0197)*IT_0201 + IT_0197*conj(IT_0201) + conj(IT_0199)*IT_0202 +
       conj(IT_0202)*(IT_0199 + IT_0205) + IT_0205*IT_0231 + IT_0232*(IT_0199 +
       2*IT_0205 + IT_0233))*IT_0234 + IT_0235*IT_0246 + (-6)*IT_0081*IT_0254 + 
      (-6)*IT_0255*IT_0257;
    return create_ccomplex_return(IT_0258);
}

