#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_sb_2_to_Z_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_sb_2_to_Z_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*V_tb*IT_0000
      *IT_0002*conj(U_sb_01)*U_st_00;
    const ccomplex_t IT_0004 = (-0.166666666666667)*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = cpow(IT_0005, -1);
    const ccomplex_t IT_0009 = cpow(IT_0006, -1);
    const ccomplex_t IT_0010 = pow(m_b, 2);
    const ccomplex_t IT_0011 = pow(m_t, 2);
    const ccomplex_t IT_0012 = cpow(IT_0005, 2);
    const ccomplex_t IT_0013 = cpow(IT_0006, 2);
    const ccomplex_t IT_0014 = m_b*mu_h*IT_0006*IT_0007*IT_0008*conj(U_sb_11)
      *U_st_00 + -conj(U_sb_01)*((IT_0007*(IT_0010 + -IT_0011) + -m_W*(IT_0012 +
       -IT_0013))*U_st_00 + m_t*mu_h*IT_0005*IT_0007*IT_0009*U_st_10);
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*V_tb*((-0.707106781186548)
      *e_em*IT_0014*IT_0016 + A_b*IT_0005*conj(U_sb_11)*U_st_00 + -A_t*IT_0006
      *conj(U_sb_01)*U_st_10);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*m_W*e_em*IT_0002*IT_0015;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = pow(m_W, 2);
    const ccomplex_t IT_0023 = pow(m_sb_2, 2);
    const ccomplex_t IT_0024 = pow(m_st_1, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + IT_0022 + -IT_0023 + -IT_0024 
      + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0021*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0023 + -IT_0024 + -reg_prop, -1);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0016*conj(U_sb_01)*U_st_00;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0016;
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0023 + 0.5*IT_0024;
    const ccomplex_t IT_0035 = s_12 + IT_0034;
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = -IT_0030;
    const ccomplex_t IT_0038 = IT_0032*IT_0037;
    const ccomplex_t IT_0039 = s_12 + IT_0023;
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = 2*IT_0032;
    const ccomplex_t IT_0042 = IT_0031*IT_0041;
    const ccomplex_t IT_0043 = s_13 + s_23;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = IT_0037*IT_0041;
    const ccomplex_t IT_0046 = s_23*IT_0045;
    const ccomplex_t IT_0047 = (-0.5)*IT_0040 + (-0.5)*IT_0044 + 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0036 + IT_0047;
    const ccomplex_t IT_0049 = IT_0029*IT_0048;
    const ccomplex_t IT_0050 = IT_0004 + -IT_0027 + (-2)*IT_0049;
    const ccomplex_t IT_0051 = pow(m_W, -2);
    const ccomplex_t IT_0052 = pow(m_Z, -2);
    const ccomplex_t IT_0053 = pow(s_34, 2);
    const ccomplex_t IT_0054 = IT_0051*IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (-2)*IT_0032;
    const ccomplex_t IT_0056 = IT_0031*IT_0055;
    const ccomplex_t IT_0057 = IT_0029*(IT_0033 + IT_0056);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = pow(s_13, 2);
    const ccomplex_t IT_0060 = IT_0052*IT_0059;
    const ccomplex_t IT_0061 = pow(s_14, 2);
    const ccomplex_t IT_0062 = IT_0051*IT_0061;
    const ccomplex_t IT_0063 = (IT_0024 + -IT_0060)*(IT_0024 + -IT_0062);
    const ccomplex_t IT_0064 = 3*IT_0058;
    const ccomplex_t IT_0065 = s_13*s_14*s_34*IT_0051*IT_0052;
    const ccomplex_t IT_0066 = -IT_0062;
    const ccomplex_t IT_0067 = IT_0024 + -IT_0060 + IT_0065 + IT_0066;
    const ccomplex_t IT_0068 = 3*IT_0050;
    const ccomplex_t IT_0069 = 3*conj(IT_0050);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0002*IT_0015 
      + (-3)*IT_0001*IT_0016)*U_st_00*conj(U_st_00) + IT_0002*IT_0015*U_st_10
      *conj(U_st_10));
    const ccomplex_t IT_0071 = 1.33333333333333*IT_0070;
    const ccomplex_t IT_0072 = IT_0037*IT_0071;
    const ccomplex_t IT_0073 = pow(m_Z, 2);
    const ccomplex_t IT_0074 = cpow((-2)*s_13 + IT_0073 + (0 + _Complex_I*1)
      *m_st_1*Gamma_t1 + reg_prop, -1);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0072*IT_0075;
    const ccomplex_t IT_0077 = IT_0033 + IT_0056;
    const ccomplex_t IT_0078 = IT_0037*IT_0055;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = IT_0077 + IT_0079;
    const ccomplex_t IT_0081 = IT_0029*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0016*conj(U_sb_01)*U_st_01;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0002*IT_0015 
      + (-3)*IT_0001*IT_0016)*U_st_00*conj(U_st_01) + IT_0002*IT_0015*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0085 = 1.33333333333333*IT_0084;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_13 + IT_0024 + IT_0073 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0086*IT_0088;
    const ccomplex_t IT_0090 = -IT_0076 + -IT_0081 + -IT_0089;
    const ccomplex_t IT_0091 = s_13*s_24*s_34*IT_0051*IT_0052;
    const ccomplex_t IT_0092 = s_14*s_24*IT_0051;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = s_13*s_23*IT_0052;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = s_12 + IT_0093 + IT_0095;
    const ccomplex_t IT_0097 = IT_0091 + IT_0096;
    const ccomplex_t IT_0098 = pow(s_24, 2);
    const ccomplex_t IT_0099 = IT_0051*IT_0098;
    const ccomplex_t IT_0100 = (IT_0024 + -IT_0060)*(IT_0023 + -IT_0099);
    const ccomplex_t IT_0101 = 3*IT_0090;
    const ccomplex_t IT_0102 = (IT_0024 + -IT_0060)*(s_12 + -IT_0092);
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = (-3)*IT_0058;
    const ccomplex_t IT_0105 = IT_0029*(IT_0033 + -IT_0038 + IT_0056);
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = s_14*s_34*IT_0051;
    const ccomplex_t IT_0108 = (s_12 + -IT_0094)*(s_13 + -IT_0107);
    const ccomplex_t IT_0109 = -IT_0032;
    const ccomplex_t IT_0110 = IT_0031*IT_0109;
    const ccomplex_t IT_0111 = IT_0029*IT_0110;
    const ccomplex_t IT_0112 = 3*IT_0111;
    const ccomplex_t IT_0113 = s_14*s_23*s_34*IT_0051*IT_0052;
    const ccomplex_t IT_0114 = IT_0096 + IT_0113;
    const ccomplex_t IT_0115 = (s_12 + -IT_0092)*(s_12 + -IT_0094);
    const ccomplex_t IT_0116 = (IT_0024 + -IT_0062)*(s_12 + -IT_0094);
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = (-3)*conj(IT_0058);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0016*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0002
      *IT_0015 + 3*IT_0001*IT_0016)*U_sb_00*conj(U_sb_01) + IT_0002*IT_0015
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0122 = (-0.666666666666667)*IT_0121;
    const ccomplex_t IT_0123 = IT_0120*IT_0122;
    const ccomplex_t IT_0124 = cpow((-2)*s_23 + IT_0023 + IT_0073 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = IT_0123*IT_0125;
    const ccomplex_t IT_0127 = -IT_0038 + -IT_0078;
    const ccomplex_t IT_0128 = IT_0077 + IT_0127;
    const ccomplex_t IT_0129 = IT_0029*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0002
      *IT_0015 + 3*IT_0001*IT_0016)*U_sb_01*conj(U_sb_01) + IT_0002*IT_0015
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0131 = (-0.666666666666667)*IT_0130;
    const ccomplex_t IT_0132 = IT_0037*IT_0131;
    const ccomplex_t IT_0133 = cpow((-2)*s_23 + IT_0073 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0132*IT_0134;
    const ccomplex_t IT_0136 = IT_0126 + -IT_0129 + IT_0135;
    const ccomplex_t IT_0137 = s_24*s_34*IT_0051;
    const ccomplex_t IT_0138 = (s_12 + -IT_0094)*(s_23 + -IT_0137);
    const ccomplex_t IT_0139 = s_23*s_24*s_34*IT_0051*IT_0052;
    const ccomplex_t IT_0140 = -IT_0099;
    const ccomplex_t IT_0141 = pow(s_23, 2);
    const ccomplex_t IT_0142 = IT_0052*IT_0141;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = IT_0023 + IT_0139 + IT_0140 + IT_0143;
    const ccomplex_t IT_0145 = (s_12 + -IT_0092)*(IT_0023 + -IT_0142);
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = (-3)*IT_0106;
    const ccomplex_t IT_0148 = (s_12 + -IT_0094)*(IT_0023 + -IT_0099);
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = (-3)*IT_0090;
    const ccomplex_t IT_0151 = (IT_0023 + -IT_0099)*(IT_0023 + -IT_0142);
    const ccomplex_t IT_0152 = IT_0037*IT_0109;
    const ccomplex_t IT_0153 = -IT_0152;
    const ccomplex_t IT_0154 = IT_0110 + IT_0153;
    const ccomplex_t IT_0155 = IT_0029*IT_0154;
    const ccomplex_t IT_0156 = -IT_0126 + -IT_0135 + IT_0155;
    const ccomplex_t IT_0157 = (s_23 + -IT_0137)*(IT_0023 + -IT_0142);
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = 3*conj(IT_0111);
    const ccomplex_t IT_0160 = (-3)*conj(IT_0106);
    const ccomplex_t IT_0161 = (-3)*conj(IT_0090);
    const ccomplex_t IT_0162 = 3*conj(IT_0058);
    const ccomplex_t IT_0163 = 0.333333333333333*IT_0115;
    const ccomplex_t IT_0164 = 3*IT_0136*(0.333333333333333*IT_0069*IT_0144 +
       conj(IT_0136)*IT_0151 + -conj(IT_0156)*IT_0158 + 0.333333333333333
      *IT_0138*IT_0159 + 0.333333333333333*IT_0146*IT_0160 + 0.333333333333333
      *IT_0149*IT_0161 + IT_0162*IT_0163);
    const ccomplex_t IT_0165 = (IT_0024 + -IT_0062)*(IT_0023 + -IT_0142);
    const ccomplex_t IT_0166 = 3*conj(IT_0090);
    const ccomplex_t IT_0167 = s_34*IT_0051;
    const ccomplex_t IT_0168 = s_13*s_34*IT_0052;
    const ccomplex_t IT_0169 = -IT_0168;
    const ccomplex_t IT_0170 = s_14 + IT_0169;
    const ccomplex_t IT_0171 = IT_0167*IT_0170;
    const ccomplex_t IT_0172 = (IT_0024 + -IT_0060)*(s_13 + -IT_0107);
    const ccomplex_t IT_0173 = -IT_0172;
    const ccomplex_t IT_0174 = (IT_0024 + -IT_0060)*(s_23 + -IT_0137);
    const ccomplex_t IT_0175 = -IT_0174;
    const ccomplex_t IT_0176 = (-0.333333333333333)*IT_0175;
    const ccomplex_t IT_0177 = IT_0051*IT_0053;
    const ccomplex_t IT_0178 = (IT_0024 + -IT_0060)*(IT_0073 + -IT_0177);
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = (-3)*conj(IT_0050);
    const ccomplex_t IT_0181 = s_23*s_34*IT_0052;
    const ccomplex_t IT_0182 = -IT_0181;
    const ccomplex_t IT_0183 = s_24 + IT_0182;
    const ccomplex_t IT_0184 = IT_0167*IT_0183;
    const ccomplex_t IT_0185 = (-0.333333333333333)*IT_0138;
    const ccomplex_t IT_0186 = (s_12 + -IT_0094)*(IT_0073 + -IT_0177);
    const ccomplex_t IT_0187 = (-0.333333333333333)*IT_0186;
    const ccomplex_t IT_0188 = (s_13 + -IT_0107)*(IT_0023 + -IT_0142);
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = (-0.333333333333333)*IT_0189;
    const ccomplex_t IT_0191 = (IT_0023 + -IT_0142)*(IT_0073 + -IT_0177);
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = 6*IT_0050*conj(IT_0050)*(1 + 0.5*IT_0054) +
       conj(IT_0058)*(IT_0063*IT_0064 + IT_0067*IT_0068) + IT_0058*IT_0067
      *IT_0069 + conj(IT_0090)*(IT_0068*IT_0097 + IT_0100*IT_0101 + IT_0103
      *IT_0104) + conj(IT_0106)*(IT_0108*IT_0112 + IT_0068*IT_0114 + IT_0101
      *IT_0115 + IT_0104*IT_0117) + IT_0090*(IT_0069*IT_0097 + IT_0103*IT_0118) 
      + conj(IT_0136)*(IT_0064*IT_0115 + IT_0112*IT_0138 + IT_0068*IT_0144 +
       IT_0146*IT_0147 + IT_0149*IT_0150) + IT_0164 + 3*IT_0106*
      (0.333333333333333*IT_0069*IT_0114 + 0.333333333333333*IT_0117*IT_0118 +
       0.333333333333333*IT_0108*IT_0159 + conj(IT_0106)*IT_0165 + IT_0163
      *IT_0166) + (-3)*conj(IT_0111)*(IT_0050*IT_0171 + (-0.333333333333333)
      *IT_0104*IT_0173 + IT_0150*IT_0176) + (-3)*IT_0111*((-0.333333333333333)
      *IT_0118*IT_0173 + IT_0161*IT_0176 + conj(IT_0111)*IT_0179 + (
      -0.333333333333333)*IT_0171*IT_0180) + (-3)*conj(IT_0156)*((
      -0.333333333333333)*IT_0064*IT_0108 + IT_0050*IT_0184 + IT_0101*IT_0185 +
       IT_0112*IT_0187 + IT_0147*IT_0190) + (-3)*IT_0156*(conj(IT_0136)*IT_0158 
      + (-0.333333333333333)*IT_0108*IT_0162 + (-0.333333333333333)*IT_0180
      *IT_0184 + IT_0166*IT_0185 + IT_0159*IT_0187 + IT_0160*IT_0190 + conj
      (IT_0156)*IT_0192);
    return create_ccomplex_return(IT_0193);
}

