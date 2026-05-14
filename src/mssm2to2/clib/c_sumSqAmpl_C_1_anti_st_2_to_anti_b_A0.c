#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_st_2_to_anti_b_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_st_2_to_anti_b_A0(
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
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
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
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = V_u1*U_Wm1;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = U_d1*V_Wp1;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (IT_0002 + IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*e_em*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, -1);
    const ccomplex_t IT_0012 = cpow(IT_0000, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0007*IT_0011*IT_0012*conj(U_st_01);
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = pow(m_b, 2);
    const ccomplex_t IT_0017 = pow(m_st_2, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = V_u1*U_Wm2;
    const ccomplex_t IT_0022 = IT_0000*IT_0021;
    const ccomplex_t IT_0023 = U_d2*V_Wp1;
    const ccomplex_t IT_0024 = IT_0003*IT_0023;
    const ccomplex_t IT_0025 = IT_0007*(IT_0022 + IT_0024);
    const ccomplex_t IT_0026 = 1.4142135623731*e_em*IT_0025;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0007*IT_0011*IT_0012*conj(U_st_01);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0020 + IT_0033;
    const ccomplex_t IT_0035 = pow(m_A0, 2);
    const ccomplex_t IT_0036 = m_b*m_C_1*IT_0035;
    const ccomplex_t IT_0037 = m_b*e_em*IT_0003*IT_0007*IT_0011*IT_0012;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0014*IT_0038;
    const ccomplex_t IT_0040 = pow(m_C_1, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + IT_0016 + -IT_0017 + -IT_0040 
      + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = 6*conj(IT_0043);
    const ccomplex_t IT_0045 = IT_0036*IT_0044;
    const ccomplex_t IT_0046 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0047 = IT_0000*IT_0046;
    const ccomplex_t IT_0048 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0049 = IT_0003*IT_0048;
    const ccomplex_t IT_0050 = IT_0007*(IT_0047 + IT_0049);
    const ccomplex_t IT_0051 = 1.4142135623731*e_em*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = IT_0031*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = cpow(IT_0003, -1);
    const ccomplex_t IT_0056 = IT_0007*IT_0055;
    const ccomplex_t IT_0057 = m_t*V_tb*V_u2*e_em*IT_0011*conj(U_st_11);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = m_C_1*IT_0061;
    const ccomplex_t IT_0063 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0064 = IT_0007*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = m_C_1*IT_0066;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0062 + IT_0068;
    const ccomplex_t IT_0070 = IT_0054*IT_0069;
    const ccomplex_t IT_0071 = m_t*V_tb*V_u1*e_em*IT_0011*conj(U_st_11);
    const ccomplex_t IT_0072 = IT_0056*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = m_C_1*IT_0075;
    const ccomplex_t IT_0077 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0078 = IT_0007*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = m_C_1*IT_0080;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = IT_0076 + IT_0082;
    const ccomplex_t IT_0084 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0085 = IT_0000*IT_0084;
    const ccomplex_t IT_0086 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0087 = IT_0003*IT_0086;
    const ccomplex_t IT_0088 = IT_0007*(IT_0085 + IT_0087);
    const ccomplex_t IT_0089 = 1.4142135623731*e_em*IT_0088;
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = IT_0018*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = IT_0083*IT_0092;
    const ccomplex_t IT_0094 = -IT_0070 + -IT_0093;
    const ccomplex_t IT_0095 = m_b*s_14;
    const ccomplex_t IT_0096 = 6*IT_0095;
    const ccomplex_t IT_0097 = conj(IT_0094)*IT_0096;
    const ccomplex_t IT_0098 = m_C_1*IT_0034;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = s_34*m_C_1;
    const ccomplex_t IT_0101 = 6*IT_0100;
    const ccomplex_t IT_0102 = conj(IT_0099)*IT_0101;
    const ccomplex_t IT_0103 = m_t*e_em*mu_h*IT_0011*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0104 = IT_0007*IT_0103;
    const ccomplex_t IT_0105 = m_t*e_em*mu_h*IT_0011*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0106 = IT_0007*IT_0105;
    const ccomplex_t IT_0107 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0108 = IT_0000*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0111 = IT_0000*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (-0.5)*IT_0104 + 0.5*IT_0106 + (-0.5)*IT_0109 +
       0.5*IT_0112;
    const ccomplex_t IT_0114 = -IT_0080;
    const ccomplex_t IT_0115 = IT_0075 + IT_0114;
    const ccomplex_t IT_0116 = IT_0113*IT_0115;
    const ccomplex_t IT_0117 = cpow((-2)*s_13 + IT_0016 + IT_0040 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0121 = IT_0010*IT_0115;
    const ccomplex_t IT_0122 = IT_0018*IT_0120*IT_0121;
    const ccomplex_t IT_0123 = cpow((-2)*s_13 + IT_0016 + IT_0040 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0124 = m_t*e_em*mu_h*IT_0011*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0125 = IT_0007*IT_0124;
    const ccomplex_t IT_0126 = m_t*e_em*mu_h*IT_0011*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0127 = IT_0007*IT_0126;
    const ccomplex_t IT_0128 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0129 = IT_0000*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0132 = IT_0000*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = (-0.5)*IT_0125 + 0.5*IT_0127 + (-0.5)*IT_0130 +
       0.5*IT_0133;
    const ccomplex_t IT_0135 = m_t*V_tb*V_u1*e_em*IT_0011*conj(U_st_10);
    const ccomplex_t IT_0136 = IT_0056*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = (-0.5)*IT_0138;
    const ccomplex_t IT_0140 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0141 = IT_0007*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = IT_0139 + IT_0144;
    const ccomplex_t IT_0146 = IT_0134*IT_0145;
    const ccomplex_t IT_0147 = IT_0123*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0150 = -IT_0066;
    const ccomplex_t IT_0151 = IT_0061 + IT_0150;
    const ccomplex_t IT_0152 = IT_0027*IT_0151;
    const ccomplex_t IT_0153 = IT_0031*IT_0149*IT_0152;
    const ccomplex_t IT_0154 = -IT_0119 + IT_0122 + -IT_0148 + IT_0153;
    const ccomplex_t IT_0155 = (-6)*IT_0095;
    const ccomplex_t IT_0156 = conj(IT_0154)*IT_0155;
    const ccomplex_t IT_0157 = IT_0014*IT_0113;
    const ccomplex_t IT_0158 = IT_0117*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0007*IT_0011*IT_0012*conj(U_st_00);
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = IT_0134*IT_0161;
    const ccomplex_t IT_0163 = IT_0123*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0029*IT_0052;
    const ccomplex_t IT_0166 = IT_0031*IT_0149*IT_0165;
    const ccomplex_t IT_0167 = IT_0014*IT_0090;
    const ccomplex_t IT_0168 = IT_0018*IT_0120*IT_0167;
    const ccomplex_t IT_0169 = -IT_0159 + -IT_0164 + IT_0166 + IT_0168;
    const ccomplex_t IT_0170 = (-6)*IT_0100;
    const ccomplex_t IT_0171 = conj(IT_0169)*IT_0170;
    const ccomplex_t IT_0172 = IT_0092*IT_0115;
    const ccomplex_t IT_0173 = IT_0054*IT_0151;
    const ccomplex_t IT_0174 = IT_0172 + IT_0173;
    const ccomplex_t IT_0175 = 6*IT_0036;
    const ccomplex_t IT_0176 = conj(IT_0174)*IT_0175;
    const ccomplex_t IT_0177 = IT_0045 + IT_0097 + IT_0102 + IT_0156 + IT_0171
       + IT_0176;
    const ccomplex_t IT_0178 = 6*IT_0043;
    const ccomplex_t IT_0179 = s_13*IT_0035;
    const ccomplex_t IT_0180 = s_14*s_34;
    const ccomplex_t IT_0181 = (-2)*IT_0180;
    const ccomplex_t IT_0182 = IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = (-6)*IT_0043;
    const ccomplex_t IT_0184 = 6*s_13;
    const ccomplex_t IT_0185 = m_b*m_C_1;
    const ccomplex_t IT_0186 = 6*IT_0185;
    const ccomplex_t IT_0187 = IT_0038*IT_0115;
    const ccomplex_t IT_0188 = IT_0041*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = conj(IT_0034) + conj(IT_0190);
    const ccomplex_t IT_0192 = (-6)*s_13;
    const ccomplex_t IT_0193 = (-6)*IT_0185;
    const ccomplex_t IT_0194 = (-6)*conj(IT_0043);
    const ccomplex_t IT_0195 = (-6)*IT_0182;
    const ccomplex_t IT_0196 = IT_0036*IT_0178;
    const ccomplex_t IT_0197 = IT_0099*IT_0101;
    const ccomplex_t IT_0198 = IT_0154*IT_0155;
    const ccomplex_t IT_0199 = IT_0034*IT_0195;
    const ccomplex_t IT_0200 = IT_0034*IT_0177 + (IT_0095*conj(IT_0099) + conj
      (IT_0094)*IT_0100)*IT_0178 + (IT_0100*conj(IT_0154) + IT_0095*conj(IT_0169
      ) + (conj(IT_0043) + conj(IT_0174))*IT_0182)*IT_0183 + IT_0099*(IT_0044
      *IT_0095 + IT_0096*conj(IT_0174) + conj(IT_0099)*IT_0184 + conj(IT_0094)
      *IT_0186 + conj(IT_0169)*IT_0192 + conj(IT_0154)*IT_0193) + IT_0094*
      (IT_0044*IT_0100 + IT_0101*conj(IT_0174) + conj(IT_0094)*IT_0184 + conj
      (IT_0099)*IT_0186 + IT_0096*IT_0191 + conj(IT_0154)*IT_0192 + conj(IT_0169
      )*IT_0193) + IT_0169*(IT_0155*conj(IT_0174) + conj(IT_0169)*IT_0184 + conj
      (IT_0154)*IT_0186 + IT_0170*IT_0191 + conj(IT_0099)*IT_0192 + conj(IT_0094
      )*IT_0193 + IT_0095*IT_0194) + IT_0154*(IT_0170*conj(IT_0174) + conj
      (IT_0154)*IT_0184 + conj(IT_0169)*IT_0186 + conj(IT_0094)*IT_0192 + conj
      (IT_0099)*IT_0193 + IT_0100*IT_0194) + IT_0174*(IT_0096*conj(IT_0099) +
       conj(IT_0094)*IT_0101 + IT_0155*conj(IT_0169) + conj(IT_0154)*IT_0170 +
       IT_0175*IT_0191 + IT_0182*IT_0194 + conj(IT_0174)*IT_0195) + IT_0190*
      (IT_0177 + IT_0191*IT_0195) + IT_0191*(IT_0196 + IT_0197 + IT_0198 +
       IT_0199);
    return create_ccomplex_return(IT_0200);
}

