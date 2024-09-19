#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_00*conj(U_st_01) + IT_0001*IT_0002*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0005 = 1.33333333333333*IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0007 = 0.333333333333333*IT_0006;
    const ccomplex_t IT_0008 = pow(m_st_1, 2);
    const ccomplex_t IT_0009 = pow(m_st_2, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0008 + -IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = IT_0005*IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = sin(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0003*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = m_t*V_tb*V_u1*e_em*IT_0018*conj(U_st_11);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0025 = IT_0003*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = IT_0023 + IT_0028;
    const ccomplex_t IT_0030 = V_tb*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0031 = IT_0003*IT_0030;
    const ccomplex_t IT_0032 = m_t*V_tb*conj(V_u1)*e_em*IT_0018*U_st_10;
    const ccomplex_t IT_0033 = IT_0017*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + (-0.5)*IT_0034);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = pow(m_b, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0009 + IT_0037 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0029*IT_0040;
    const ccomplex_t IT_0042 = m_t*V_tb*V_u2*e_em*IT_0018*conj(U_st_11);
    const ccomplex_t IT_0043 = IT_0017*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0048 = IT_0003*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0046 + IT_0051;
    const ccomplex_t IT_0053 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0054 = IT_0003*IT_0053;
    const ccomplex_t IT_0055 = m_t*V_tb*conj(V_u2)*e_em*IT_0018*U_st_10;
    const ccomplex_t IT_0056 = IT_0017*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0054 + (-0.5)*IT_0057);
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_23 + IT_0009 + IT_0037 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0052*IT_0062;
    const ccomplex_t IT_0064 = IT_0041 + IT_0063;
    const ccomplex_t IT_0065 = IT_0001*IT_0002;
    const ccomplex_t IT_0066 = e_em*IT_0065;
    const ccomplex_t IT_0067 = IT_0000*IT_0003;
    const ccomplex_t IT_0068 = e_em*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0066 + 3*IT_0068);
    const ccomplex_t IT_0070 = (-0.166666666666667)*IT_0069;
    const ccomplex_t IT_0071 = IT_0010*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0005*IT_0072;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = cos(beta);
    const ccomplex_t IT_0076 = cpow(IT_0075, -1);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0003*IT_0018*IT_0076*U_st_00;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0003*IT_0018*IT_0076*conj(U_st_01);
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0078*IT_0080;
    const ccomplex_t IT_0082 = IT_0038*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0003*IT_0018*IT_0076*U_st_00;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0003*IT_0018*IT_0076*conj(U_st_01);
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = IT_0060*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0083 + IT_0090;
    const ccomplex_t IT_0092 = IT_0008*IT_0037;
    const ccomplex_t IT_0093 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0094 = m_b*IT_0093;
    const ccomplex_t IT_0095 = m_b*IT_0005;
    const ccomplex_t IT_0096 = IT_0012*(IT_0094 + -IT_0095);
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = IT_0072*(IT_0094 + -IT_0095);
    const ccomplex_t IT_0099 = -IT_0072*(IT_0094 + -IT_0095);
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = cos(alpha);
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*m_b*e_em*IT_0003*IT_0018
      *IT_0076*IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = sin(alpha);
    const ccomplex_t IT_0106 = m_t*mu_h*IT_0018;
    const ccomplex_t IT_0107 = IT_0075*IT_0102;
    const ccomplex_t IT_0108 = -IT_0015*IT_0105 + IT_0107;
    const ccomplex_t IT_0109 = cpow(IT_0000, -2);
    const ccomplex_t IT_0110 = pow(m_t, 2);
    const ccomplex_t IT_0111 = (-2)*U_st_00*((0.5*m_W*IT_0108*(IT_0003 + -1./3
      *IT_0002*IT_0109) + IT_0003*IT_0016*IT_0018*IT_0105*IT_0110)*conj(U_st_01)
       + -1./2*IT_0003*IT_0016*IT_0102*IT_0106*conj(U_st_11)) + U_st_10*(IT_0003
      *IT_0016*IT_0102*IT_0106*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0108*IT_0109 + 1.5*IT_0003*IT_0016*IT_0018*IT_0105*IT_0110)*conj
      (U_st_11));
    const ccomplex_t IT_0112 = (0 + _Complex_I*-1)*(e_em*IT_0111 + (
      -1.4142135623731)*A_t*IT_0105*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0104*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_b*e_em*IT_0003*IT_0018
      *IT_0076*IT_0105;
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0075*IT_0105;
    const ccomplex_t IT_0121 = IT_0015*IT_0102;
    const ccomplex_t IT_0122 = IT_0120 + IT_0121;
    const ccomplex_t IT_0123 = 2*U_st_00*((IT_0003*IT_0016*IT_0018*IT_0102
      *IT_0110 + (-0.5)*m_W*(IT_0003 + -1./3*IT_0002*IT_0109)*IT_0122)*conj
      (U_st_01) + 1./2*IT_0003*IT_0016*IT_0105*IT_0106*conj(U_st_11)) + U_st_10*
      (IT_0003*IT_0016*IT_0105*IT_0106*conj(U_st_01) + (-1.33333333333333)*((
      -1.5)*IT_0003*IT_0016*IT_0018*IT_0102*IT_0110 + m_W*IT_0002*IT_0109
      *IT_0122)*conj(U_st_11));
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(e_em*IT_0123 +
       1.4142135623731*A_t*IT_0102*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = IT_0119*IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0008 + -IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0128 = IT_0126*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = -IT_0117 + -IT_0129;
    const ccomplex_t IT_0131 = 18*s_34;
    const ccomplex_t IT_0132 = (-18)*IT_0037;
    const ccomplex_t IT_0133 = IT_0131 + IT_0132;
    const ccomplex_t IT_0134 = s_34*IT_0008;
    const ccomplex_t IT_0135 = (-6)*IT_0134;
    const ccomplex_t IT_0136 = s_13*s_14;
    const ccomplex_t IT_0137 = 12*IT_0136;
    const ccomplex_t IT_0138 = IT_0135 + IT_0137;
    const ccomplex_t IT_0139 = (-18)*IT_0134;
    const ccomplex_t IT_0140 = 36*IT_0136;
    const ccomplex_t IT_0141 = IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = IT_0036*IT_0080;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0144 = IT_0038*IT_0142*IT_0143;
    const ccomplex_t IT_0145 = IT_0059*IT_0087;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0147 = IT_0060*IT_0145*IT_0146;
    const ccomplex_t IT_0148 = m_b*IT_0083;
    const ccomplex_t IT_0149 = m_b*IT_0090;
    const ccomplex_t IT_0150 = -IT_0144 + -IT_0147 + -IT_0148 + -IT_0149;
    const ccomplex_t IT_0151 = (-6)*IT_0101;
    const ccomplex_t IT_0152 = (-0.166666666666667)*IT_0151;
    const ccomplex_t IT_0153 = IT_0029*IT_0078;
    const ccomplex_t IT_0154 = IT_0038*IT_0143*IT_0153;
    const ccomplex_t IT_0155 = IT_0052*IT_0085;
    const ccomplex_t IT_0156 = IT_0060*IT_0146*IT_0155;
    const ccomplex_t IT_0157 = -IT_0154 + -IT_0156;
    const ccomplex_t IT_0158 = (-18)*IT_0157;
    const ccomplex_t IT_0159 = (-0.166666666666667)*IT_0158;
    const ccomplex_t IT_0160 = m_b*IT_0023;
    const ccomplex_t IT_0161 = m_b*IT_0027;
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = IT_0160 + IT_0162;
    const ccomplex_t IT_0164 = IT_0040*IT_0163;
    const ccomplex_t IT_0165 = m_b*IT_0046;
    const ccomplex_t IT_0166 = m_b*IT_0050;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = IT_0062*IT_0168;
    const ccomplex_t IT_0170 = -IT_0164 + -IT_0169;
    const ccomplex_t IT_0171 = (-18)*IT_0170;
    const ccomplex_t IT_0172 = (-0.166666666666667)*IT_0171;
    const ccomplex_t IT_0173 = IT_0096 + IT_0152 + IT_0159 + IT_0172;
    const ccomplex_t IT_0174 = conj(IT_0097) + conj(IT_0098);
    const ccomplex_t IT_0175 = (-18)*IT_0096;
    const ccomplex_t IT_0176 = (-18)*IT_0101;
    const ccomplex_t IT_0177 = (-6)*IT_0170;
    const ccomplex_t IT_0178 = (-6)*IT_0157;
    const ccomplex_t IT_0179 = IT_0175 + IT_0176 + IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = IT_0097 + IT_0098;
    const ccomplex_t IT_0181 = (-18)*conj(IT_0096);
    const ccomplex_t IT_0182 = (-18)*conj(IT_0101);
    const ccomplex_t IT_0183 = (-6)*conj(IT_0170);
    const ccomplex_t IT_0184 = (-6)*conj(IT_0157);
    const ccomplex_t IT_0185 = IT_0181 + IT_0182 + IT_0183 + IT_0184;
    const ccomplex_t IT_0186 = (-18)*conj(IT_0170);
    const ccomplex_t IT_0187 = (-18)*conj(IT_0157);
    const ccomplex_t IT_0188 = (-6)*conj(IT_0096);
    const ccomplex_t IT_0189 = (-6)*conj(IT_0101);
    const ccomplex_t IT_0190 = IT_0186 + IT_0187 + IT_0188 + IT_0189;
    const ccomplex_t IT_0191 = m_b*s_14;
    const ccomplex_t IT_0192 = 3*IT_0150;
    const ccomplex_t IT_0193 = IT_0180 + IT_0192;
    const ccomplex_t IT_0194 = 3*conj(IT_0150);
    const ccomplex_t IT_0195 = IT_0174 + IT_0194;
    const ccomplex_t IT_0196 = conj(IT_0091)*IT_0173 + IT_0074*(conj(IT_0150) 
      + 3*IT_0174) + (-0.166666666666667)*conj(IT_0014)*IT_0179 + conj(IT_0074)*
      (IT_0150 + 3*IT_0180) + (-0.166666666666667)*IT_0014*IT_0185 + (
      -0.166666666666667)*IT_0091*IT_0190 + conj(IT_0064)*IT_0193 + IT_0064
      *IT_0195;
    const ccomplex_t IT_0197 = conj(IT_0096) + conj(IT_0101);
    const ccomplex_t IT_0198 = IT_0157 + IT_0170;
    const ccomplex_t IT_0199 = conj(IT_0157) + conj(IT_0170);
    const ccomplex_t IT_0200 = IT_0096 + IT_0101;
    const ccomplex_t IT_0201 = 6*s_34;
    const ccomplex_t IT_0202 = m_b*s_13;
    const ccomplex_t IT_0203 = 6*IT_0150;
    const ccomplex_t IT_0204 = 6*conj(IT_0150);
    const ccomplex_t IT_0205 = conj(IT_0064)*(IT_0096 + IT_0101 + 3*IT_0157 +
       3*IT_0170) + 0.166666666666667*conj(IT_0074)*(18*IT_0096 + 18*IT_0101 + 6
      *IT_0157 + 6*IT_0170) + IT_0064*(conj(IT_0096) + conj(IT_0101) + 3*conj
      (IT_0157) + 3*conj(IT_0170)) + 0.166666666666667*IT_0074*(18*conj(IT_0096)
       + 18*conj(IT_0101) + 6*conj(IT_0157) + 6*conj(IT_0170)) + conj(IT_0091)
      *IT_0193 + IT_0091*IT_0195 + 0.166666666666667*conj(IT_0014)*(18*IT_0097 +
       18*IT_0098 + IT_0203) + 0.166666666666667*IT_0014*(18*conj(IT_0097) + 18
      *conj(IT_0098) + IT_0204);
    const ccomplex_t IT_0206 = -IT_0037;
    const ccomplex_t IT_0207 = 6*IT_0130;
    const ccomplex_t IT_0208 = 6*conj(IT_0130);
    const ccomplex_t IT_0209 = IT_0097*IT_0131;
    const ccomplex_t IT_0210 = IT_0101*IT_0131;
    const ccomplex_t IT_0211 = -s_14;
    const ccomplex_t IT_0212 = s_13 + IT_0211;
    const ccomplex_t IT_0213 = m_b*IT_0212;
    const ccomplex_t IT_0214 = 18*IT_0130;
    const ccomplex_t IT_0215 = 18*conj(IT_0130);
    const ccomplex_t IT_0216 = pow(m_Z, 2);
    const ccomplex_t IT_0217 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + IT_0216 
      + -reg_prop, -1);
    const ccomplex_t IT_0218 = m_b*e_em*IT_0003*IT_0018;
    const ccomplex_t IT_0219 = 0.5*IT_0218;
    const ccomplex_t IT_0220 = m_t*e_em*mu_h*IT_0018*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0221 = IT_0003*IT_0016*IT_0075;
    const ccomplex_t IT_0222 = IT_0220*IT_0221;
    const ccomplex_t IT_0223 = m_t*e_em*mu_h*IT_0018*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0224 = IT_0221*IT_0223;
    const ccomplex_t IT_0225 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0226 = IT_0015*IT_0225;
    const ccomplex_t IT_0227 = 1.4142135623731*IT_0226;
    const ccomplex_t IT_0228 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0229 = IT_0015*IT_0228;
    const ccomplex_t IT_0230 = 1.4142135623731*IT_0229;
    const ccomplex_t IT_0231 = 0.5*IT_0222 + (-0.5)*IT_0224 + (-0.5)*IT_0227 +
       0.5*IT_0230;
    const ccomplex_t IT_0232 = IT_0219*IT_0231;
    const ccomplex_t IT_0233 = IT_0217*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0236 = m_b*e_em*IT_0003*IT_0015*IT_0018*IT_0076;
    const ccomplex_t IT_0237 = (-0.5)*IT_0236;
    const ccomplex_t IT_0238 = IT_0003*IT_0220;
    const ccomplex_t IT_0239 = IT_0003*IT_0223;
    const ccomplex_t IT_0240 = IT_0075*IT_0225;
    const ccomplex_t IT_0241 = 1.4142135623731*IT_0240;
    const ccomplex_t IT_0242 = IT_0075*IT_0228;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = (-0.5)*IT_0238 + 0.5*IT_0239 + (-0.5)*IT_0241 +
       0.5*IT_0243;
    const ccomplex_t IT_0245 = IT_0237*IT_0244;
    const ccomplex_t IT_0246 = IT_0235*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = -IT_0234 + -IT_0247;
    const ccomplex_t IT_0249 = s_13 + s_14;
    const ccomplex_t IT_0250 = m_b*IT_0249;
    const ccomplex_t IT_0251 = (-18)*s_34;
    const ccomplex_t IT_0252 = IT_0132 + IT_0251;
    const ccomplex_t IT_0253 = 18*IT_0037;
    const ccomplex_t IT_0254 = IT_0131 + IT_0253;
    const ccomplex_t IT_0255 = s_34 + IT_0037;
    const ccomplex_t IT_0256 = 36*s_34;
    const ccomplex_t IT_0257 = (18*conj(IT_0014) + (-6)*conj(IT_0064) + (-18)
      *conj(IT_0074) + 6*conj(IT_0091))*IT_0250 + IT_0197*IT_0252 + IT_0174
      *IT_0254 + (IT_0183 + IT_0184 + IT_0204)*IT_0255 + conj(IT_0248)*(36
      *IT_0037 + IT_0256);
    const ccomplex_t IT_0258 = (-6)*(conj(IT_0014)*(IT_0064 + 3*IT_0074) +
       IT_0014*(conj(IT_0064) + 3*conj(IT_0074)) + 3*(conj(IT_0064) +
       0.333333333333333*conj(IT_0074))*IT_0091 + 3*(IT_0064 + 0.333333333333333
      *IT_0074)*conj(IT_0091))*IT_0092 + ((conj(IT_0096) + conj(IT_0097) + conj
      (IT_0098) + conj(IT_0101))*IT_0130 + (IT_0096 + IT_0097 + IT_0098 +
       IT_0101)*conj(IT_0130))*IT_0133 + (conj(IT_0064)*IT_0074 + IT_0064*conj
      (IT_0074) + conj(IT_0014)*IT_0091 + IT_0014*conj(IT_0091))*IT_0138 + 
      (IT_0014*conj(IT_0014) + IT_0064*conj(IT_0064) + IT_0074*conj(IT_0074) +
       IT_0091*conj(IT_0091))*IT_0141 + (-6)*IT_0037*(conj(IT_0150)*IT_0173 + (
      -0.166666666666667)*IT_0174*IT_0179 + (-0.166666666666667)*IT_0180*IT_0185
       + (-0.166666666666667)*IT_0150*IT_0190) + (-6)*IT_0191*IT_0196 + IT_0131*
      (IT_0150*conj(IT_0150) + IT_0098*IT_0174 + IT_0096*IT_0197 + IT_0198
      *IT_0199) + (IT_0150*IT_0174 + conj(IT_0150)*IT_0180 + IT_0197*IT_0198 +
       IT_0199*IT_0200)*IT_0201 + 6*IT_0202*IT_0205 + (s_34 + IT_0206)*((conj
      (IT_0150) + IT_0199)*IT_0207 + (IT_0150 + IT_0198)*IT_0208) + IT_0174
      *IT_0209 + IT_0197*IT_0210 + IT_0213*((conj(IT_0064) + conj(IT_0091))
      *IT_0207 + (IT_0064 + IT_0091)*IT_0208 + (conj(IT_0014) + conj(IT_0074))
      *IT_0214 + (IT_0014 + IT_0074)*IT_0215) + conj(IT_0248)*((18*IT_0014 + (-6
      )*IT_0064 + (-18)*IT_0074 + 6*IT_0091)*IT_0250 + IT_0200*IT_0252 + IT_0180
      *IT_0254 + (IT_0177 + IT_0178 + IT_0203)*IT_0255) + IT_0130*conj(IT_0130)*
      ((-36)*IT_0037 + IT_0256) + IT_0248*IT_0257;
    return create_ccomplex_return(IT_0258);
}

