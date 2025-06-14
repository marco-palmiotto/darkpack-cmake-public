#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_st_1_to_anti_b_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_st_1_to_anti_b_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
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
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = m_b*m_C_2;
    const ccomplex_t IT_0001 = pow(m_b, 2);
    const ccomplex_t IT_0002 = pow(m_C_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0004*IT_0006*IT_0008*conj(U_st_00);
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = m_t*mu_h*IT_0004;
    const ccomplex_t IT_0013 = cos(alpha);
    const ccomplex_t IT_0014 = sin(beta);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = IT_0005*IT_0013;
    const ccomplex_t IT_0017 = -IT_0011*IT_0014 + IT_0016;
    const ccomplex_t IT_0018 = cos(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -2);
    const ccomplex_t IT_0020 = pow(m_t, 2);
    const ccomplex_t IT_0021 = (-2)*U_st_00*((0.5*m_W*IT_0017*(IT_0008 + -1./3
      *IT_0007*IT_0019) + IT_0004*IT_0008*IT_0011*IT_0015*IT_0020)*conj(U_st_00)
       + -1./2*IT_0008*IT_0012*IT_0013*IT_0015*conj(U_st_10)) + U_st_10*(IT_0008
      *IT_0012*IT_0013*IT_0015*conj(U_st_00) + (-1.33333333333333)*(m_W*IT_0007
      *IT_0017*IT_0019 + 1.5*IT_0004*IT_0008*IT_0011*IT_0015*IT_0020)*conj
      (U_st_10));
    const ccomplex_t IT_0022 = (0 + _Complex_I*-1)*(e_em*IT_0021 + (
      -1.4142135623731)*A_t*IT_0011*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0010*IT_0023;
    const ccomplex_t IT_0025 = IT_0003*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0028 = (-2)*conj(U_st_00)*((0.5*m_W*IT_0017*(IT_0008 +
       -1./3*IT_0007*IT_0019) + IT_0004*IT_0008*IT_0011*IT_0015*IT_0020)*U_st_01
       + -1./2*IT_0008*IT_0012*IT_0013*IT_0015*U_st_11) + conj(U_st_10)*(IT_0008
      *IT_0012*IT_0013*IT_0015*U_st_01 + (-1.33333333333333)*(m_W*IT_0007
      *IT_0017*IT_0019 + 1.5*IT_0004*IT_0008*IT_0011*IT_0015*IT_0020)*U_st_11);
    const ccomplex_t IT_0029 = (0 + _Complex_I*-1)*(e_em*IT_0028 + (
      -1.4142135623731)*A_t*IT_0011*(U_st_01*conj(U_st_10) + conj(U_st_00)
      *U_st_11));
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0004*IT_0006*IT_0008*conj(U_st_01);
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0037 = IT_0013*IT_0036;
    const ccomplex_t IT_0038 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0039 = IT_0011*IT_0038;
    const ccomplex_t IT_0040 = IT_0008*(IT_0037 + IT_0039);
    const ccomplex_t IT_0041 = 1.4142135623731*e_em*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = IT_0010*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0046 = pow(m_st_1, 2);
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0001 + IT_0046 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0044*IT_0045*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*m_b*e_em*IT_0004*IT_0006
      *IT_0008*IT_0013;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0010*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_12 + IT_0001 + -IT_0002 + -IT_0046 
      + -reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = m_b*IT_0054;
    const ccomplex_t IT_0056 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0057 = IT_0013*IT_0056;
    const ccomplex_t IT_0058 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0059 = IT_0011*IT_0058;
    const ccomplex_t IT_0060 = IT_0008*(IT_0057 + IT_0059);
    const ccomplex_t IT_0061 = 1.4142135623731*e_em*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0004*IT_0006*IT_0008*conj(U_st_00);
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0068 = cpow((-2)*s_23 + IT_0001 + IT_0046 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0066*IT_0067*IT_0068;
    const ccomplex_t IT_0070 = -IT_0026 + -IT_0035 + IT_0048 + 2*IT_0055 +
       IT_0069;
    const ccomplex_t IT_0071 = IT_0043*IT_0047;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0008*IT_0015;
    const ccomplex_t IT_0074 = m_t*V_tb*V_u2*e_em*IT_0004*conj(U_st_10);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = m_C_2*IT_0078;
    const ccomplex_t IT_0080 = V_tb*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0081 = IT_0008*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = m_C_2*IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0079 + IT_0085;
    const ccomplex_t IT_0087 = IT_0072*IT_0086;
    const ccomplex_t IT_0088 = IT_0063*IT_0068;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = m_t*V_tb*V_u1*e_em*IT_0004*conj(U_st_10);
    const ccomplex_t IT_0091 = IT_0073*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = m_C_2*IT_0094;
    const ccomplex_t IT_0096 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0097 = IT_0008*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = m_C_2*IT_0099;
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = IT_0095 + IT_0101;
    const ccomplex_t IT_0103 = IT_0089*IT_0102;
    const ccomplex_t IT_0104 = -IT_0087 + -IT_0103;
    const ccomplex_t IT_0105 = U_d1*V_Wp2;
    const ccomplex_t IT_0106 = IT_0013*IT_0105;
    const ccomplex_t IT_0107 = V_u2*U_Wm1;
    const ccomplex_t IT_0108 = IT_0011*IT_0107;
    const ccomplex_t IT_0109 = IT_0008*(IT_0106 + IT_0108);
    const ccomplex_t IT_0110 = 1.4142135623731*e_em*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = -IT_0099;
    const ccomplex_t IT_0114 = IT_0094 + IT_0113;
    const ccomplex_t IT_0115 = IT_0112*IT_0114;
    const ccomplex_t IT_0116 = IT_0067*IT_0068*IT_0115;
    const ccomplex_t IT_0117 = IT_0050*IT_0052;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = m_b*IT_0078;
    const ccomplex_t IT_0120 = m_b*IT_0083;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = IT_0119 + IT_0121;
    const ccomplex_t IT_0123 = IT_0118*IT_0122;
    const ccomplex_t IT_0124 = -IT_0083;
    const ccomplex_t IT_0125 = IT_0078 + IT_0124;
    const ccomplex_t IT_0126 = IT_0023*IT_0125;
    const ccomplex_t IT_0127 = IT_0003*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = U_d2*V_Wp2;
    const ccomplex_t IT_0130 = IT_0013*IT_0129;
    const ccomplex_t IT_0131 = V_u2*U_Wm2;
    const ccomplex_t IT_0132 = IT_0011*IT_0131;
    const ccomplex_t IT_0133 = IT_0008*(IT_0130 + IT_0132);
    const ccomplex_t IT_0134 = 1.4142135623731*e_em*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = 0.5*IT_0135;
    const ccomplex_t IT_0137 = IT_0125*IT_0136;
    const ccomplex_t IT_0138 = IT_0045*IT_0047*IT_0137;
    const ccomplex_t IT_0139 = m_t*V_tb*V_u2*e_em*IT_0004*conj(U_st_11);
    const ccomplex_t IT_0140 = IT_0073*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0145 = IT_0008*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = -IT_0146;
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = IT_0143 + IT_0148;
    const ccomplex_t IT_0150 = IT_0030*IT_0149;
    const ccomplex_t IT_0151 = IT_0027*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = IT_0116 + 2*IT_0123 + -IT_0128 + IT_0138 + 
      -IT_0152;
    const ccomplex_t IT_0154 = IT_0010*IT_0136;
    const ccomplex_t IT_0155 = IT_0047*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = IT_0065*IT_0112;
    const ccomplex_t IT_0158 = IT_0068*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = IT_0156 + IT_0159;
    const ccomplex_t IT_0161 = m_C_2*IT_0160;
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = 6*IT_0104;
    const ccomplex_t IT_0164 = (-6)*IT_0153;
    const ccomplex_t IT_0165 = IT_0163 + IT_0164;
    const ccomplex_t IT_0166 = 6*conj(IT_0104);
    const ccomplex_t IT_0167 = (-6)*conj(IT_0153);
    const ccomplex_t IT_0168 = IT_0166 + IT_0167;
    const ccomplex_t IT_0169 = m_b*s_14;
    const ccomplex_t IT_0170 = -IT_0054;
    const ccomplex_t IT_0171 = IT_0072*IT_0125;
    const ccomplex_t IT_0172 = IT_0089*IT_0114;
    const ccomplex_t IT_0173 = IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = conj(IT_0170) + conj(IT_0173);
    const ccomplex_t IT_0175 = 6*IT_0162;
    const ccomplex_t IT_0176 = (-6)*IT_0070;
    const ccomplex_t IT_0177 = IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = IT_0170 + IT_0173;
    const ccomplex_t IT_0179 = 6*conj(IT_0162);
    const ccomplex_t IT_0180 = (-6)*conj(IT_0070);
    const ccomplex_t IT_0181 = IT_0179 + IT_0180;
    const ccomplex_t IT_0182 = (IT_0078 + -IT_0083)*IT_0118;
    const ccomplex_t IT_0183 = -IT_0182;
    const ccomplex_t IT_0184 = conj(IT_0160) + conj(IT_0183);
    const ccomplex_t IT_0185 = IT_0160 + IT_0183;
    const ccomplex_t IT_0186 = s_34*m_C_2;
    const ccomplex_t IT_0187 = 6*s_13;
    const ccomplex_t IT_0188 = pow(m_H0, 2);
    const ccomplex_t IT_0189 = m_b*m_C_2*IT_0188;
    const ccomplex_t IT_0190 = 6*IT_0160;
    const ccomplex_t IT_0191 = 6*conj(IT_0160);
    const ccomplex_t IT_0192 = (-6)*s_13;
    const ccomplex_t IT_0193 = s_13*IT_0188;
    const ccomplex_t IT_0194 = s_14*s_34;
    const ccomplex_t IT_0195 = (-2)*IT_0194;
    const ccomplex_t IT_0196 = IT_0193 + IT_0195;
    const ccomplex_t IT_0197 = (-6)*IT_0160;
    const ccomplex_t IT_0198 = (-6)*IT_0173;
    const ccomplex_t IT_0199 = (-6)*IT_0000*(conj(IT_0070)*(IT_0104 + -IT_0153
      ) + IT_0070*(conj(IT_0104) + -conj(IT_0153)) + (-0.166666666666667)*conj
      (IT_0162)*IT_0165 + (-0.166666666666667)*IT_0162*IT_0168) + IT_0169*
      (IT_0174*IT_0177 + IT_0178*IT_0181 + IT_0165*IT_0184 + IT_0168*IT_0185) + 
      (IT_0165*IT_0174 + IT_0168*IT_0178 + IT_0177*IT_0184 + IT_0181*IT_0185)
      *IT_0186 + (IT_0070*conj(IT_0070) + IT_0104*conj(IT_0104) + IT_0153*conj
      (IT_0153) + IT_0162*conj(IT_0162))*IT_0187 + 6*IT_0189*(IT_0174*(IT_0183 +
       0.166666666666667*IT_0190) + IT_0178*(conj(IT_0183) + 0.166666666666667
      *IT_0191)) + (conj(IT_0104)*IT_0153 + IT_0104*conj(IT_0153) + conj(IT_0070
      )*IT_0162 + IT_0070*conj(IT_0162))*IT_0192 + (-6)*IT_0196*(IT_0184*
      (IT_0183 + (-0.166666666666667)*IT_0197) + IT_0174*(IT_0170 + (
      -0.166666666666667)*IT_0198));
    return create_ccomplex_return(IT_0199);
}

