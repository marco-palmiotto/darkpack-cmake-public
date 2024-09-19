#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_anti_sb_2_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_anti_sb_2_to_Z_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = U_sb_01*conj(U_sb_01);
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0001, -4);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = IT_0005*IT_0007;
    const ccomplex_t IT_0009 = IT_0003*IT_0008;
    const ccomplex_t IT_0010 = U_sb_11*conj(U_sb_11);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0006, -2);
    const ccomplex_t IT_0013 = IT_0003*IT_0012;
    const ccomplex_t IT_0014 = tan(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, 2);
    const ccomplex_t IT_0016 = cpow(1 + IT_0015, -1);
    const ccomplex_t IT_0017 = (IT_0004 + 0.166666666666667*IT_0009 +
       0.666666666666667*IT_0011 + 1.5*IT_0013)*IT_0016;
    const ccomplex_t IT_0018 = IT_0000*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = 0.333333333333333*IT_0019;
    const ccomplex_t IT_0021 = m_W*e_em;
    const ccomplex_t IT_0022 = cos(alpha);
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = IT_0005*IT_0006*IT_0022*IT_0023;
    const ccomplex_t IT_0025 = cos(beta);
    const ccomplex_t IT_0026 = sin(alpha);
    const ccomplex_t IT_0027 = IT_0005*IT_0006*IT_0025*IT_0026;
    const ccomplex_t IT_0028 = cpow(IT_0006, -3);
    const ccomplex_t IT_0029 = IT_0022*IT_0023*IT_0028;
    const ccomplex_t IT_0030 = IT_0025*IT_0026*IT_0028;
    const ccomplex_t IT_0031 = cpow(IT_0006, -1);
    const ccomplex_t IT_0032 = IT_0002*IT_0022*IT_0023*IT_0031;
    const ccomplex_t IT_0033 = IT_0002*IT_0025*IT_0026*IT_0031;
    const ccomplex_t IT_0034 = cpow(IT_0002 + IT_0012, -1);
    const ccomplex_t IT_0035 = (IT_0024 + -IT_0027 + IT_0029 + -IT_0030 + 2
      *IT_0032 + (-2)*IT_0033)*IT_0034;
    const ccomplex_t IT_0036 = IT_0021*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow(V_tb, 2);
    const ccomplex_t IT_0039 = pow(m_W, -1);
    const ccomplex_t IT_0040 = m_b*mu_h*IT_0039;
    const ccomplex_t IT_0041 = cpow(IT_0025, -1);
    const ccomplex_t IT_0042 = IT_0025*IT_0026;
    const ccomplex_t IT_0043 = IT_0022*IT_0023;
    const ccomplex_t IT_0044 = IT_0042 + IT_0043;
    const ccomplex_t IT_0045 = pow(m_b, 2);
    const ccomplex_t IT_0046 = 2*U_sb_01*(((-0.166666666666667)*m_W*(IT_0002
      *IT_0006 + 3*IT_0031)*IT_0044 + IT_0026*IT_0031*IT_0039*IT_0041*IT_0045)
      *conj(U_sb_01) + 1./2*IT_0022*IT_0031*IT_0040*IT_0041*conj(U_sb_11)) +
       U_sb_11*(IT_0022*IT_0031*IT_0040*IT_0041*conj(U_sb_01) + (
      -0.666666666666667)*(m_W*IT_0002*IT_0006*IT_0044 + (-3)*IT_0026*IT_0031
      *IT_0039*IT_0041*IT_0045)*conj(U_sb_11));
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(e_em*IT_0046 +
       1.4142135623731*A_b*IT_0026*IT_0038*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0037*IT_0048;
    const ccomplex_t IT_0050 = pow(m_sb_2, 2);
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0050 + -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0005*IT_0006*IT_0022*IT_0025;
    const ccomplex_t IT_0055 = IT_0005*IT_0006*IT_0023*IT_0026;
    const ccomplex_t IT_0056 = IT_0022*IT_0025*IT_0028;
    const ccomplex_t IT_0057 = IT_0023*IT_0026*IT_0028;
    const ccomplex_t IT_0058 = IT_0002*IT_0022*IT_0025*IT_0031;
    const ccomplex_t IT_0059 = IT_0002*IT_0023*IT_0026*IT_0031;
    const ccomplex_t IT_0060 = IT_0034*(IT_0054 + IT_0055 + IT_0056 + IT_0057 
      + 2*IT_0058 + 2*IT_0059);
    const ccomplex_t IT_0061 = IT_0021*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = -IT_0023*IT_0026;
    const ccomplex_t IT_0064 = IT_0022*IT_0025;
    const ccomplex_t IT_0065 = IT_0063 + IT_0064;
    const ccomplex_t IT_0066 = (-2)*U_sb_01*((IT_0022*IT_0031*IT_0039*IT_0041
      *IT_0045 + (-0.166666666666667)*m_W*(IT_0002*IT_0006 + 3*IT_0031)*IT_0065)
      *conj(U_sb_01) + -1./2*IT_0026*IT_0031*IT_0040*IT_0041*conj(U_sb_11)) +
       U_sb_11*(IT_0026*IT_0031*IT_0040*IT_0041*conj(U_sb_01) + (-2)*(IT_0022
      *IT_0031*IT_0039*IT_0041*IT_0045 + -1./3*m_W*IT_0002*IT_0006*IT_0065)*conj
      (U_sb_11));
    const ccomplex_t IT_0067 = (0 + _Complex_I*-1)*(e_em*IT_0066 + (
      -1.4142135623731)*A_b*IT_0022*IT_0038*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0062*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_12 + (-2)*IT_0050 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0020 + -IT_0053 + -IT_0072;
    const ccomplex_t IT_0074 = pow(m_Z, -4);
    const ccomplex_t IT_0075 = pow(s_34, 2);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = cpow(IT_0001, -1);
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0001
      *IT_0031 + 1./3*IT_0006*IT_0077)*U_sb_01*conj(U_sb_01) + IT_0006*IT_0077
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0079 = (-0.333333333333333)*IT_0078;
    const ccomplex_t IT_0080 = (-0.666666666666667)*IT_0078;
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = pow(m_Z, 2);
    const ccomplex_t IT_0083 = cpow((-2)*s_13 + IT_0082 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = IT_0081*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0001
      *IT_0031 + 1./3*IT_0006*IT_0077)*U_sb_00*conj(U_sb_01) + IT_0006*IT_0077
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0087 = (-0.333333333333333)*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0001
      *IT_0031 + 1./3*IT_0006*IT_0077)*conj(U_sb_00)*U_sb_01 + IT_0006*IT_0077
      *conj(U_sb_10)*U_sb_11);
    const ccomplex_t IT_0089 = (-0.666666666666667)*IT_0088;
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_13 + IT_0050 + IT_0082 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = IT_0085 + IT_0093;
    const ccomplex_t IT_0095 = pow(m_Z, -2);
    const ccomplex_t IT_0096 = pow(s_13, 2);
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = -IT_0050;
    const ccomplex_t IT_0099 = pow(s_14, 2);
    const ccomplex_t IT_0100 = IT_0095*IT_0099;
    const ccomplex_t IT_0101 = (IT_0097 + IT_0098)*(IT_0098 + IT_0100);
    const ccomplex_t IT_0102 = 3*IT_0094;
    const ccomplex_t IT_0103 = s_13*s_14*s_34*IT_0074;
    const ccomplex_t IT_0104 = IT_0097 + IT_0098 + IT_0100 + -IT_0103;
    const ccomplex_t IT_0105 = (-3)*IT_0073;
    const ccomplex_t IT_0106 = (-3)*conj(IT_0073);
    const ccomplex_t IT_0107 = cpow((-2)*s_23 + IT_0082 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = cpow(IT_0080, 2);
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = IT_0081 + IT_0110;
    const ccomplex_t IT_0112 = IT_0108*IT_0111;
    const ccomplex_t IT_0113 = cpow((-2)*s_23 + IT_0050 + IT_0082 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = (-0.666666666666667)*IT_0086;
    const ccomplex_t IT_0116 = IT_0089*IT_0115;
    const ccomplex_t IT_0117 = (-0.333333333333333)*IT_0088;
    const ccomplex_t IT_0118 = IT_0115*IT_0117;
    const ccomplex_t IT_0119 = -IT_0118;
    const ccomplex_t IT_0120 = IT_0116 + IT_0119;
    const ccomplex_t IT_0121 = IT_0114*IT_0120;
    const ccomplex_t IT_0122 = IT_0112 + -IT_0121;
    const ccomplex_t IT_0123 = s_13*s_23*IT_0095;
    const ccomplex_t IT_0124 = s_14*s_34*IT_0095;
    const ccomplex_t IT_0125 = (s_12 + -IT_0123)*(s_13 + -IT_0124);
    const ccomplex_t IT_0126 = 3*conj(IT_0094);
    const ccomplex_t IT_0127 = s_34*IT_0095;
    const ccomplex_t IT_0128 = s_23*s_34*IT_0095;
    const ccomplex_t IT_0129 = -IT_0128;
    const ccomplex_t IT_0130 = s_24 + IT_0129;
    const ccomplex_t IT_0131 = IT_0127*IT_0130;
    const ccomplex_t IT_0132 = pow(s_23, 2);
    const ccomplex_t IT_0133 = IT_0095*IT_0132;
    const ccomplex_t IT_0134 = (s_13 + -IT_0124)*(IT_0050 + -IT_0133);
    const ccomplex_t IT_0135 = -IT_0134;
    const ccomplex_t IT_0136 = IT_0081*IT_0108;
    const ccomplex_t IT_0137 = IT_0114*IT_0118;
    const ccomplex_t IT_0138 = -IT_0136 + -IT_0137;
    const ccomplex_t IT_0139 = (-3)*conj(IT_0138);
    const ccomplex_t IT_0140 = s_14*s_23*s_34*IT_0074;
    const ccomplex_t IT_0141 = s_14*s_24*IT_0095;
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = -IT_0123;
    const ccomplex_t IT_0144 = s_12 + IT_0142 + IT_0143;
    const ccomplex_t IT_0145 = IT_0140 + IT_0144;
    const ccomplex_t IT_0146 = 3*IT_0073;
    const ccomplex_t IT_0147 = -(IT_0050 + -IT_0100)*(IT_0050 + -IT_0133);
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = 3*IT_0138;
    const ccomplex_t IT_0150 = -(IT_0050 + -IT_0100)*(s_12 + -IT_0123);
    const ccomplex_t IT_0151 = (-3)*IT_0094;
    const ccomplex_t IT_0152 = 3*conj(IT_0073);
    const ccomplex_t IT_0153 = (-3)*conj(IT_0094);
    const ccomplex_t IT_0154 = IT_0084*IT_0111;
    const ccomplex_t IT_0155 = -IT_0116;
    const ccomplex_t IT_0156 = IT_0090 + IT_0155;
    const ccomplex_t IT_0157 = IT_0092*IT_0156;
    const ccomplex_t IT_0158 = IT_0154 + IT_0157;
    const ccomplex_t IT_0159 = s_13*s_24*s_34*IT_0074;
    const ccomplex_t IT_0160 = IT_0144 + IT_0159;
    const ccomplex_t IT_0161 = (s_12 + -IT_0123)*(s_12 + -IT_0141);
    const ccomplex_t IT_0162 = 3*conj(IT_0138);
    const ccomplex_t IT_0163 = -(IT_0050 + -IT_0097)*(s_12 + -IT_0141);
    const ccomplex_t IT_0164 = -IT_0112 + IT_0121;
    const ccomplex_t IT_0165 = s_23*s_24*s_34*IT_0074;
    const ccomplex_t IT_0166 = pow(s_24, 2);
    const ccomplex_t IT_0167 = IT_0095*IT_0166;
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = -IT_0133;
    const ccomplex_t IT_0170 = IT_0050 + IT_0165 + IT_0168 + IT_0169;
    const ccomplex_t IT_0171 = (IT_0050 + -IT_0133)*(s_12 + -IT_0141);
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = (s_12 + -IT_0123)*(IT_0050 + -IT_0167);
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = (-3)*conj(IT_0158);
    const ccomplex_t IT_0176 = -IT_0085 + -IT_0093;
    const ccomplex_t IT_0177 = s_13*s_34*IT_0095;
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = s_14 + IT_0178;
    const ccomplex_t IT_0180 = IT_0127*IT_0179;
    const ccomplex_t IT_0181 = -(IT_0050 + -IT_0097)*(s_13 + -IT_0124);
    const ccomplex_t IT_0182 = s_24*s_34*IT_0095;
    const ccomplex_t IT_0183 = -(IT_0050 + -IT_0097)*(s_23 + -IT_0182);
    const ccomplex_t IT_0184 = (s_12 + -IT_0123)*(s_23 + -IT_0182);
    const ccomplex_t IT_0185 = -(IT_0050 + -IT_0097)*(IT_0050 + -IT_0167);
    const ccomplex_t IT_0186 = -IT_0185;
    const ccomplex_t IT_0187 = 3*IT_0158;
    const ccomplex_t IT_0188 = 0.333333333333333*IT_0161;
    const ccomplex_t IT_0189 = (IT_0050 + -IT_0133)*(IT_0050 + -IT_0167);
    const ccomplex_t IT_0190 = 3*IT_0164;
    const ccomplex_t IT_0191 = (-3)*IT_0138;
    const ccomplex_t IT_0192 = (-3)*IT_0158;
    const ccomplex_t IT_0193 = (IT_0050 + -IT_0133)*(s_23 + -IT_0182);
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = -IT_0122;
    const ccomplex_t IT_0196 = IT_0075*IT_0095;
    const ccomplex_t IT_0197 = (s_12 + -IT_0123)*(IT_0082 + -IT_0196);
    const ccomplex_t IT_0198 = (IT_0050 + -IT_0133)*(IT_0082 + -IT_0196);
    const ccomplex_t IT_0199 = -IT_0198;
    const ccomplex_t IT_0200 = (-0.333333333333333)*IT_0184;
    const ccomplex_t IT_0201 = -(IT_0050 + -IT_0097)*(IT_0082 + -IT_0196);
    const ccomplex_t IT_0202 = 6*IT_0073*conj(IT_0073)*(1 + 0.5*IT_0076) +
       conj(IT_0094)*(IT_0101*IT_0102 + IT_0104*IT_0105) + IT_0094*IT_0104
      *IT_0106 + IT_0122*(IT_0125*IT_0126 + IT_0106*IT_0131 + IT_0135*IT_0139) +
       conj(IT_0138)*(IT_0145*IT_0146 + IT_0148*IT_0149 + IT_0150*IT_0151) +
       IT_0138*(IT_0145*IT_0152 + IT_0150*IT_0153) + IT_0158*(IT_0152*IT_0160 +
       IT_0161*IT_0162 + IT_0153*IT_0163) + IT_0164*(IT_0126*IT_0161 + IT_0152
      *IT_0170 + IT_0139*IT_0172 + IT_0174*IT_0175) + IT_0176*(IT_0125*IT_0162 +
       IT_0106*IT_0180 + IT_0153*IT_0181 + IT_0175*IT_0183) + 3*conj(IT_0158)*
      (0.333333333333333*IT_0146*IT_0160 + 0.333333333333333*IT_0151*IT_0163 +
       IT_0122*IT_0184 + 0.333333333333333*IT_0186*IT_0187 + IT_0149*IT_0188) +
       3*conj(IT_0164)*(0.333333333333333*IT_0146*IT_0170 + IT_0176*IT_0184 +
       IT_0102*IT_0188 + 0.333333333333333*IT_0189*IT_0190 + 0.333333333333333
      *IT_0172*IT_0191 + 0.333333333333333*IT_0174*IT_0192 + IT_0194*IT_0195) + 
      (-3)*conj(IT_0122)*((-0.333333333333333)*IT_0102*IT_0125 + (
      -0.333333333333333)*IT_0105*IT_0131 + (-0.333333333333333)*IT_0135*IT_0191
       + IT_0164*IT_0194 + -IT_0176*IT_0197 + IT_0122*IT_0199 + IT_0187*IT_0200)
       + (-3)*conj(IT_0176)*((-0.333333333333333)*IT_0125*IT_0149 + (
      -0.333333333333333)*IT_0105*IT_0180 + (-0.333333333333333)*IT_0151*IT_0181
       + (-0.333333333333333)*IT_0183*IT_0192 + IT_0195*IT_0197 + IT_0190
      *IT_0200 + IT_0176*IT_0201);
    return create_ccomplex_return(IT_0202);
}

