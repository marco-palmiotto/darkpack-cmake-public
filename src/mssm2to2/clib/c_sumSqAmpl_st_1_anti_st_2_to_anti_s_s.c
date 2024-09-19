#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_anti_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
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
    const ccomplex_t V_ts = param->V_ts;
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
    const ccomplex_t IT_0015 = conj(V_ts)*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0016 = IT_0003*IT_0015;
    const ccomplex_t IT_0017 = sin(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0003*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = m_t*conj(V_ts)*conj(V_u1)*e_em*IT_0020*U_st_10;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0016 + (-0.5)*IT_0023);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = pow(m_s, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_23 + IT_0009 + IT_0026 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = m_t*V_ts*V_u1*e_em*IT_0020*conj(U_st_11);
    const ccomplex_t IT_0031 = IT_0019*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = V_ts*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0036 = IT_0003*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = IT_0034 + IT_0039;
    const ccomplex_t IT_0041 = IT_0029*IT_0040;
    const ccomplex_t IT_0042 = m_t*V_ts*V_u2*e_em*IT_0020*conj(U_st_11);
    const ccomplex_t IT_0043 = IT_0019*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = V_ts*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0048 = IT_0003*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0046 + IT_0051;
    const ccomplex_t IT_0053 = conj(V_ts)*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0054 = IT_0003*IT_0053;
    const ccomplex_t IT_0055 = m_t*conj(V_ts)*conj(V_u2)*e_em*IT_0020*U_st_10;
    const ccomplex_t IT_0056 = IT_0019*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0054 + (-0.5)*IT_0057);
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_23 + IT_0009 + IT_0026 + -m_C_2*
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
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_ts)*e_em*IT_0003*IT_0020*IT_0076*U_st_00;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0003*IT_0020*IT_0076*conj(U_st_01);
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0078*IT_0080;
    const ccomplex_t IT_0082 = IT_0027*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0003*IT_0020*IT_0076*U_st_00;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0003*IT_0020*IT_0076*conj(U_st_01);
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = IT_0060*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0083 + IT_0090;
    const ccomplex_t IT_0092 = IT_0008*IT_0026;
    const ccomplex_t IT_0093 = s_34*IT_0008;
    const ccomplex_t IT_0094 = (-6)*IT_0093;
    const ccomplex_t IT_0095 = s_13*s_14;
    const ccomplex_t IT_0096 = 12*IT_0095;
    const ccomplex_t IT_0097 = IT_0094 + IT_0096;
    const ccomplex_t IT_0098 = (-18)*IT_0093;
    const ccomplex_t IT_0099 = 36*IT_0095;
    const ccomplex_t IT_0100 = IT_0098 + IT_0099;
    const ccomplex_t IT_0101 = 18*s_34;
    const ccomplex_t IT_0102 = (-18)*IT_0026;
    const ccomplex_t IT_0103 = IT_0101 + IT_0102;
    const ccomplex_t IT_0104 = sin(alpha);
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*m_s*e_em*IT_0003*IT_0020
      *IT_0076*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = cos(alpha);
    const ccomplex_t IT_0108 = m_t*mu_h*IT_0020;
    const ccomplex_t IT_0109 = IT_0075*IT_0104;
    const ccomplex_t IT_0110 = IT_0017*IT_0107;
    const ccomplex_t IT_0111 = IT_0109 + IT_0110;
    const ccomplex_t IT_0112 = cpow(IT_0000, -2);
    const ccomplex_t IT_0113 = pow(m_t, 2);
    const ccomplex_t IT_0114 = 2*U_st_00*(((-0.5)*m_W*IT_0111*(IT_0003 + -1./3
      *IT_0002*IT_0112) + IT_0003*IT_0018*IT_0020*IT_0107*IT_0113)*conj(U_st_01)
       + 1./2*IT_0003*IT_0018*IT_0104*IT_0108*conj(U_st_11)) + U_st_10*(IT_0003
      *IT_0018*IT_0104*IT_0108*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0111*IT_0112 + (-1.5)*IT_0003*IT_0018*IT_0020*IT_0107*IT_0113)*conj
      (U_st_11));
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(e_em*IT_0114 +
       1.4142135623731*A_t*IT_0107*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = IT_0106*IT_0116;
    const ccomplex_t IT_0118 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0008 + -IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0119 = IT_0117*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*m_s*e_em*IT_0003*IT_0020
      *IT_0076*IT_0107;
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = IT_0075*IT_0107;
    const ccomplex_t IT_0124 = -IT_0017*IT_0104 + IT_0123;
    const ccomplex_t IT_0125 = (-2)*U_st_00*((IT_0003*IT_0018*IT_0020*IT_0104
      *IT_0113 + 0.5*m_W*(IT_0003 + -1./3*IT_0002*IT_0112)*IT_0124)*conj(U_st_01
      ) + -1./2*IT_0003*IT_0018*IT_0107*IT_0108*conj(U_st_11)) + U_st_10*
      (IT_0003*IT_0018*IT_0107*IT_0108*conj(U_st_01) + (-1.33333333333333)*(1.5
      *IT_0003*IT_0018*IT_0020*IT_0104*IT_0113 + m_W*IT_0002*IT_0112*IT_0124)
      *conj(U_st_11));
    const ccomplex_t IT_0126 = (0 + _Complex_I*-1)*(e_em*IT_0125 + (
      -1.4142135623731)*A_t*IT_0104*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0127 = (-0.5)*IT_0126;
    const ccomplex_t IT_0128 = IT_0122*IT_0127;
    const ccomplex_t IT_0129 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0130 = IT_0128*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = -IT_0120 + -IT_0131;
    const ccomplex_t IT_0133 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0134 = m_s*IT_0133;
    const ccomplex_t IT_0135 = m_s*IT_0005;
    const ccomplex_t IT_0136 = IT_0072*(IT_0134 + -IT_0135);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = -IT_0012*(IT_0134 + -IT_0135);
    const ccomplex_t IT_0139 = -IT_0138;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = IT_0136 + IT_0140;
    const ccomplex_t IT_0142 = conj(IT_0136) + conj(IT_0140);
    const ccomplex_t IT_0143 = m_s*s_14;
    const ccomplex_t IT_0144 = IT_0025*IT_0080;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0146 = IT_0027*IT_0144*IT_0145;
    const ccomplex_t IT_0147 = IT_0059*IT_0087;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0149 = IT_0060*IT_0147*IT_0148;
    const ccomplex_t IT_0150 = m_s*IT_0083;
    const ccomplex_t IT_0151 = m_s*IT_0090;
    const ccomplex_t IT_0152 = -IT_0146 + -IT_0149 + -IT_0150 + -IT_0151;
    const ccomplex_t IT_0153 = 3*IT_0152;
    const ccomplex_t IT_0154 = IT_0141 + IT_0153;
    const ccomplex_t IT_0155 = 3*conj(IT_0152);
    const ccomplex_t IT_0156 = IT_0142 + IT_0155;
    const ccomplex_t IT_0157 = m_s*IT_0046;
    const ccomplex_t IT_0158 = m_s*IT_0050;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = IT_0157 + IT_0159;
    const ccomplex_t IT_0161 = IT_0062*IT_0160;
    const ccomplex_t IT_0162 = m_s*IT_0034;
    const ccomplex_t IT_0163 = m_s*IT_0038;
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = IT_0162 + IT_0164;
    const ccomplex_t IT_0166 = IT_0029*IT_0165;
    const ccomplex_t IT_0167 = -IT_0161 + -IT_0166;
    const ccomplex_t IT_0168 = (-18)*IT_0167;
    const ccomplex_t IT_0169 = IT_0052*IT_0085;
    const ccomplex_t IT_0170 = IT_0060*IT_0148*IT_0169;
    const ccomplex_t IT_0171 = IT_0040*IT_0078;
    const ccomplex_t IT_0172 = IT_0027*IT_0145*IT_0171;
    const ccomplex_t IT_0173 = -IT_0170 + -IT_0172;
    const ccomplex_t IT_0174 = (-18)*IT_0173;
    const ccomplex_t IT_0175 = (-6)*IT_0137;
    const ccomplex_t IT_0176 = (-6)*IT_0139;
    const ccomplex_t IT_0177 = IT_0168 + IT_0174 + IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = (-18)*conj(IT_0167);
    const ccomplex_t IT_0179 = (-18)*conj(IT_0173);
    const ccomplex_t IT_0180 = (-6)*conj(IT_0137);
    const ccomplex_t IT_0181 = (-6)*conj(IT_0139);
    const ccomplex_t IT_0182 = IT_0178 + IT_0179 + IT_0180 + IT_0181;
    const ccomplex_t IT_0183 = (-18)*IT_0137;
    const ccomplex_t IT_0184 = (-18)*IT_0139;
    const ccomplex_t IT_0185 = (-6)*IT_0167;
    const ccomplex_t IT_0186 = (-6)*IT_0173;
    const ccomplex_t IT_0187 = IT_0183 + IT_0184 + IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = (-18)*conj(IT_0137);
    const ccomplex_t IT_0189 = (-18)*conj(IT_0139);
    const ccomplex_t IT_0190 = (-6)*conj(IT_0167);
    const ccomplex_t IT_0191 = (-6)*conj(IT_0173);
    const ccomplex_t IT_0192 = IT_0188 + IT_0189 + IT_0190 + IT_0191;
    const ccomplex_t IT_0193 = 3*conj(IT_0074)*(IT_0141 + 1./3*IT_0152) + 3
      *IT_0074*(IT_0142 + 1./3*conj(IT_0152)) + conj(IT_0064)*IT_0154 + IT_0064
      *IT_0156 + (-0.166666666666667)*conj(IT_0091)*IT_0177 + (
      -0.166666666666667)*IT_0091*IT_0182 + (-0.166666666666667)*conj(IT_0014)
      *IT_0187 + (-0.166666666666667)*IT_0014*IT_0192;
    const ccomplex_t IT_0194 = IT_0167 + IT_0173;
    const ccomplex_t IT_0195 = conj(IT_0167) + conj(IT_0173);
    const ccomplex_t IT_0196 = conj(IT_0137) + conj(IT_0139);
    const ccomplex_t IT_0197 = IT_0137 + IT_0139;
    const ccomplex_t IT_0198 = 6*s_34;
    const ccomplex_t IT_0199 = -IT_0026;
    const ccomplex_t IT_0200 = 6*IT_0132;
    const ccomplex_t IT_0201 = 6*conj(IT_0132);
    const ccomplex_t IT_0202 = m_s*s_13;
    const ccomplex_t IT_0203 = 6*IT_0152;
    const ccomplex_t IT_0204 = 6*conj(IT_0152);
    const ccomplex_t IT_0205 = 6*conj(IT_0091)*IT_0154 + 6*IT_0091*IT_0156 + 6
      *conj(IT_0064)*(IT_0137 + IT_0139 + 3*IT_0167 + 3*IT_0173) + 18*conj
      (IT_0074)*(IT_0137 + IT_0139 + 1./3*IT_0167 + 1./3*IT_0173) + 6*IT_0064*
      (conj(IT_0137) + conj(IT_0139) + 3*conj(IT_0167) + 3*conj(IT_0173)) + 18
      *IT_0074*(conj(IT_0137) + conj(IT_0139) + 1./3*conj(IT_0167) + 1./3*conj
      (IT_0173)) + conj(IT_0014)*(18*IT_0136 + 18*IT_0140 + IT_0203) + IT_0014*
      (18*conj(IT_0136) + 18*conj(IT_0140) + IT_0204);
    const ccomplex_t IT_0206 = IT_0101*IT_0137;
    const ccomplex_t IT_0207 = -s_14;
    const ccomplex_t IT_0208 = s_13 + IT_0207;
    const ccomplex_t IT_0209 = m_s*IT_0208;
    const ccomplex_t IT_0210 = 18*IT_0132;
    const ccomplex_t IT_0211 = 18*conj(IT_0132);
    const ccomplex_t IT_0212 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0213 = m_s*e_em*IT_0003*IT_0017*IT_0020*IT_0076;
    const ccomplex_t IT_0214 = (-0.5)*IT_0213;
    const ccomplex_t IT_0215 = m_t*e_em*mu_h*IT_0020*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0216 = IT_0003*IT_0215;
    const ccomplex_t IT_0217 = m_t*e_em*mu_h*IT_0020*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0218 = IT_0003*IT_0217;
    const ccomplex_t IT_0219 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0220 = IT_0075*IT_0219;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0223 = IT_0075*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = (-0.5)*IT_0216 + 0.5*IT_0218 + (-0.5)*IT_0221 +
       0.5*IT_0224;
    const ccomplex_t IT_0226 = IT_0214*IT_0225;
    const ccomplex_t IT_0227 = IT_0212*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = pow(m_Z, 2);
    const ccomplex_t IT_0230 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + IT_0229 
      + -reg_prop, -1);
    const ccomplex_t IT_0231 = m_s*e_em*IT_0003*IT_0020;
    const ccomplex_t IT_0232 = 0.5*IT_0231;
    const ccomplex_t IT_0233 = IT_0003*IT_0018*IT_0075;
    const ccomplex_t IT_0234 = IT_0215*IT_0233;
    const ccomplex_t IT_0235 = IT_0217*IT_0233;
    const ccomplex_t IT_0236 = IT_0017*IT_0219;
    const ccomplex_t IT_0237 = 1.4142135623731*IT_0236;
    const ccomplex_t IT_0238 = IT_0017*IT_0222;
    const ccomplex_t IT_0239 = 1.4142135623731*IT_0238;
    const ccomplex_t IT_0240 = 0.5*IT_0234 + (-0.5)*IT_0235 + (-0.5)*IT_0237 +
       0.5*IT_0239;
    const ccomplex_t IT_0241 = IT_0232*IT_0240;
    const ccomplex_t IT_0242 = IT_0230*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = -IT_0228 + -IT_0243;
    const ccomplex_t IT_0245 = s_13 + s_14;
    const ccomplex_t IT_0246 = m_s*IT_0245;
    const ccomplex_t IT_0247 = 18*IT_0026;
    const ccomplex_t IT_0248 = IT_0101 + IT_0247;
    const ccomplex_t IT_0249 = (-18)*s_34;
    const ccomplex_t IT_0250 = IT_0102 + IT_0249;
    const ccomplex_t IT_0251 = s_34 + IT_0026;
    const ccomplex_t IT_0252 = 36*s_34;
    const ccomplex_t IT_0253 = (18*conj(IT_0014) + (-6)*conj(IT_0064) + (-18)
      *conj(IT_0074) + 6*conj(IT_0091))*IT_0246 + IT_0142*IT_0248 + IT_0196
      *IT_0250 + (IT_0190 + IT_0191 + IT_0204)*IT_0251 + conj(IT_0244)*(36
      *IT_0026 + IT_0252);
    const ccomplex_t IT_0254 = (-6)*(conj(IT_0014)*(IT_0064 + 3*IT_0074) +
       IT_0014*(conj(IT_0064) + 3*conj(IT_0074)) + 3*(conj(IT_0064) +
       0.333333333333333*conj(IT_0074))*IT_0091 + 3*(IT_0064 + 0.333333333333333
      *IT_0074)*conj(IT_0091))*IT_0092 + (conj(IT_0064)*IT_0074 + IT_0064*conj
      (IT_0074) + conj(IT_0014)*IT_0091 + IT_0014*conj(IT_0091))*IT_0097 + 
      (IT_0014*conj(IT_0014) + IT_0064*conj(IT_0064) + IT_0074*conj(IT_0074) +
       IT_0091*conj(IT_0091))*IT_0100 + IT_0103*(conj(IT_0132)*(IT_0137 +
       IT_0139 + IT_0141) + IT_0132*(conj(IT_0137) + conj(IT_0139) + IT_0142)) +
       IT_0026*(conj(IT_0152)*IT_0177 + IT_0152*IT_0182 + IT_0142*IT_0187 +
       IT_0141*IT_0192) + (-6)*IT_0143*IT_0193 + IT_0101*(IT_0141*IT_0142 +
       IT_0152*conj(IT_0152) + IT_0194*IT_0195 + IT_0139*IT_0196) + (IT_0142
      *IT_0152 + IT_0141*conj(IT_0152) + IT_0194*IT_0196 + IT_0195*IT_0197)
      *IT_0198 + (s_34 + IT_0199)*((conj(IT_0152) + IT_0195)*IT_0200 + (IT_0152 
      + IT_0194)*IT_0201) + IT_0202*IT_0205 + IT_0196*IT_0206 + IT_0209*((conj
      (IT_0064) + conj(IT_0091))*IT_0200 + (IT_0064 + IT_0091)*IT_0201 + (conj
      (IT_0014) + conj(IT_0074))*IT_0210 + (IT_0014 + IT_0074)*IT_0211) + conj
      (IT_0244)*((18*IT_0014 + (-6)*IT_0064 + (-18)*IT_0074 + 6*IT_0091)*IT_0246
       + IT_0141*IT_0248 + IT_0197*IT_0250 + (IT_0185 + IT_0186 + IT_0203)
      *IT_0251) + IT_0132*conj(IT_0132)*((-36)*IT_0026 + IT_0252) + IT_0244
      *IT_0253;
    return create_ccomplex_return(IT_0254);
}

