#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
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
    const creal_t m_sb_1 = param->m_sb_1;
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
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_00*conj(U_sb_00) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0005 = (-0.333333333333333)*IT_0004;
    const ccomplex_t IT_0006 = m_t*IT_0005;
    const ccomplex_t IT_0007 = (-0.666666666666667)*IT_0004;
    const ccomplex_t IT_0008 = m_t*IT_0007;
    const ccomplex_t IT_0009 = pow(m_sb_1, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0001*IT_0002;
    const ccomplex_t IT_0012 = e_em*IT_0011;
    const ccomplex_t IT_0013 = IT_0000*IT_0003;
    const ccomplex_t IT_0014 = e_em*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0012 + (-3)*IT_0014);
    const ccomplex_t IT_0016 = (-0.166666666666667)*IT_0015;
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0020 = (-0.666666666666667)*IT_0019;
    const ccomplex_t IT_0021 = IT_0010*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (IT_0006 + -IT_0008)*(IT_0018 + -IT_0022);
    const ccomplex_t IT_0024 = m_t*s_13;
    const ccomplex_t IT_0025 = IT_0007*IT_0022;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = cos(alpha);
    const ccomplex_t IT_0029 = sin(beta);
    const ccomplex_t IT_0030 = cpow(IT_0029, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*m_t*e_em*IT_0003*IT_0027
      *IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = cpow(V_tb, 2);
    const ccomplex_t IT_0034 = sin(alpha);
    const ccomplex_t IT_0035 = m_b*mu_h*IT_0027;
    const ccomplex_t IT_0036 = cos(beta);
    const ccomplex_t IT_0037 = cpow(IT_0036, -1);
    const ccomplex_t IT_0038 = IT_0034*IT_0036;
    const ccomplex_t IT_0039 = IT_0028*IT_0029;
    const ccomplex_t IT_0040 = IT_0038 + IT_0039;
    const ccomplex_t IT_0041 = cpow(IT_0000, -2);
    const ccomplex_t IT_0042 = pow(m_b, 2);
    const ccomplex_t IT_0043 = 2*U_sb_00*(((-0.5)*m_W*IT_0040*(IT_0003 + 1./3
      *IT_0002*IT_0041) + IT_0003*IT_0027*IT_0034*IT_0037*IT_0042)*conj(U_sb_00)
       + 1./2*IT_0003*IT_0028*IT_0035*IT_0037*conj(U_sb_10)) + U_sb_10*(IT_0003
      *IT_0028*IT_0035*IT_0037*conj(U_sb_00) + (-0.666666666666667)*(m_W*IT_0002
      *IT_0040*IT_0041 + (-3)*IT_0003*IT_0027*IT_0034*IT_0037*IT_0042)*conj
      (U_sb_10));
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(e_em*IT_0043 +
       1.4142135623731*A_b*IT_0033*IT_0034*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = IT_0032*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_12 + (-2)*IT_0009 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*m_t*e_em*IT_0003*IT_0027
      *IT_0030*IT_0034;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = -IT_0029*IT_0034;
    const ccomplex_t IT_0054 = IT_0028*IT_0036;
    const ccomplex_t IT_0055 = IT_0053 + IT_0054;
    const ccomplex_t IT_0056 = (-2)*U_sb_00*((IT_0003*IT_0027*IT_0028*IT_0037
      *IT_0042 + (-0.5)*m_W*(IT_0003 + 1./3*IT_0002*IT_0041)*IT_0055)*conj
      (U_sb_00) + -1./2*IT_0003*IT_0034*IT_0035*IT_0037*conj(U_sb_10)) + U_sb_10
      *(IT_0003*IT_0034*IT_0035*IT_0037*conj(U_sb_00) + (-2)*(IT_0003*IT_0027
      *IT_0028*IT_0037*IT_0042 + -1./3*m_W*IT_0002*IT_0041*IT_0055)*conj(U_sb_10
      ));
    const ccomplex_t IT_0057 = (0 + _Complex_I*-1)*(e_em*IT_0056 + (
      -1.4142135623731)*A_b*IT_0028*IT_0033*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0052*IT_0058;
    const ccomplex_t IT_0060 = IT_0050*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = -IT_0049 + -IT_0061;
    const ccomplex_t IT_0063 = pow(m_t, 2);
    const ccomplex_t IT_0064 = (-18)*IT_0063;
    const ccomplex_t IT_0065 = 18*s_34;
    const ccomplex_t IT_0066 = IT_0064 + IT_0065;
    const ccomplex_t IT_0067 = m_t*e_em*IT_0003*IT_0027;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = m_b*e_em*mu_h*IT_0027*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0070 = IT_0003*IT_0029*IT_0037;
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = m_b*e_em*mu_h*IT_0027*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0073 = IT_0070*IT_0072;
    const ccomplex_t IT_0074 = A_b*IT_0033*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0075 = IT_0036*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = A_b*IT_0033*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0078 = IT_0036*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0071 + 0.5*IT_0073 + 0.5*IT_0076 + (
      -0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0068*IT_0080;
    const ccomplex_t IT_0082 = pow(m_Z, 2);
    const ccomplex_t IT_0083 = cpow((-2)*s_12 + (-2)*IT_0009 + IT_0082 + 
      -reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_12 + (-2)*IT_0009 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0003*IT_0069;
    const ccomplex_t IT_0088 = IT_0003*IT_0072;
    const ccomplex_t IT_0089 = IT_0029*IT_0074;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = IT_0029*IT_0077;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (-0.5)*IT_0087 + 0.5*IT_0088 + (-0.5)*IT_0090 +
       0.5*IT_0092;
    const ccomplex_t IT_0094 = m_t*e_em*IT_0003*IT_0027*IT_0030*IT_0036;
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0093*IT_0095;
    const ccomplex_t IT_0097 = IT_0086*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = -IT_0085 + -IT_0098;
    const ccomplex_t IT_0100 = 18*IT_0063;
    const ccomplex_t IT_0101 = IT_0065 + IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u1)*e_em*IT_0003*IT_0027*IT_0030*conj(U_sb_00);
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0105 = IT_0003*IT_0104;
    const ccomplex_t IT_0106 = IT_0003*IT_0037;
    const ccomplex_t IT_0107 = m_b*conj(U_d1)*V_tb*e_em*IT_0027*U_sb_10;
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0105 + (-0.5)*IT_0109);
    const ccomplex_t IT_0111 = IT_0103*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0113 = cpow((-2)*s_23 + IT_0009 + IT_0063 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0114 = IT_0111*IT_0112*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0003*IT_0027*IT_0030*conj(U_sb_00);
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0118 = IT_0003*IT_0117;
    const ccomplex_t IT_0119 = m_b*conj(U_d2)*V_tb*e_em*IT_0027*U_sb_10;
    const ccomplex_t IT_0120 = IT_0106*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0118 + (-0.5)*IT_0121);
    const ccomplex_t IT_0123 = IT_0116*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0125 = cpow((-2)*s_23 + IT_0009 + IT_0063 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0126 = IT_0123*IT_0124*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0003*IT_0027*IT_0030*U_sb_00;
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = IT_0103*IT_0128;
    const ccomplex_t IT_0130 = IT_0113*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = m_t*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0003*IT_0027*IT_0030*U_sb_00;
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = IT_0116*IT_0134;
    const ccomplex_t IT_0136 = IT_0125*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = m_t*IT_0137;
    const ccomplex_t IT_0139 = -IT_0114 + -IT_0126 + IT_0132 + IT_0138;
    const ccomplex_t IT_0140 = 6*s_34;
    const ccomplex_t IT_0141 = -IT_0023;
    const ccomplex_t IT_0142 = 36*s_34;
    const ccomplex_t IT_0143 = V_tb*e_em*U_Wm1*conj(U_sb_00);
    const ccomplex_t IT_0144 = IT_0003*IT_0143;
    const ccomplex_t IT_0145 = m_b*U_d1*V_tb*e_em*IT_0027*conj(U_sb_10);
    const ccomplex_t IT_0146 = IT_0106*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*(IT_0144 + (-0.5)*IT_0147);
    const ccomplex_t IT_0149 = IT_0128*IT_0148;
    const ccomplex_t IT_0150 = IT_0112*IT_0113*IT_0149;
    const ccomplex_t IT_0151 = V_tb*e_em*U_Wm2*conj(U_sb_00);
    const ccomplex_t IT_0152 = IT_0003*IT_0151;
    const ccomplex_t IT_0153 = m_b*U_d2*V_tb*e_em*IT_0027*conj(U_sb_10);
    const ccomplex_t IT_0154 = IT_0106*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0152 + (-0.5)*IT_0155);
    const ccomplex_t IT_0157 = IT_0134*IT_0156;
    const ccomplex_t IT_0158 = IT_0124*IT_0125*IT_0157;
    const ccomplex_t IT_0159 = IT_0110*IT_0148;
    const ccomplex_t IT_0160 = IT_0113*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = m_t*IT_0161;
    const ccomplex_t IT_0163 = IT_0122*IT_0156;
    const ccomplex_t IT_0164 = IT_0125*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = m_t*IT_0165;
    const ccomplex_t IT_0167 = -IT_0150 + -IT_0158 + IT_0162 + IT_0166;
    const ccomplex_t IT_0168 = IT_0139 + IT_0167;
    const ccomplex_t IT_0169 = -IT_0063;
    const ccomplex_t IT_0170 = s_34 + IT_0169;
    const ccomplex_t IT_0171 = 6*IT_0170;
    const ccomplex_t IT_0172 = (-18)*conj(IT_0023);
    const ccomplex_t IT_0173 = (-6)*conj(IT_0023);
    const ccomplex_t IT_0174 = (-18)*s_34 + IT_0064;
    const ccomplex_t IT_0175 = -IT_0161 + -IT_0165;
    const ccomplex_t IT_0176 = 6*IT_0024;
    const ccomplex_t IT_0177 = IT_0007*IT_0018;
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = 18*IT_0024;
    const ccomplex_t IT_0180 = (-18)*IT_0023;
    const ccomplex_t IT_0181 = m_t*s_14;
    const ccomplex_t IT_0182 = (-18)*IT_0181;
    const ccomplex_t IT_0183 = (-6)*IT_0063;
    const ccomplex_t IT_0184 = -IT_0131 + -IT_0137;
    const ccomplex_t IT_0185 = (-6)*IT_0181;
    const ccomplex_t IT_0186 = conj(IT_0141)*(IT_0062*IT_0066 + IT_0065
      *IT_0141 + IT_0140*IT_0167 + IT_0099*IT_0174 + IT_0175*IT_0176 + IT_0178
      *IT_0179 + IT_0063*IT_0180 + IT_0026*IT_0182 + IT_0139*IT_0183 + IT_0184
      *IT_0185);
    const ccomplex_t IT_0187 = s_34 + IT_0063;
    const ccomplex_t IT_0188 = 6*IT_0187;
    const ccomplex_t IT_0189 = s_13 + s_14;
    const ccomplex_t IT_0190 = m_t*IT_0189;
    const ccomplex_t IT_0191 = 6*IT_0190;
    const ccomplex_t IT_0192 = 18*IT_0190;
    const ccomplex_t IT_0193 = (-18)*IT_0190;
    const ccomplex_t IT_0194 = (-6)*IT_0187;
    const ccomplex_t IT_0195 = (-6)*IT_0190;
    const ccomplex_t IT_0196 = conj(IT_0099)*(IT_0023*IT_0101 + IT_0099*(36
      *IT_0063 + IT_0142) + IT_0141*IT_0174 + IT_0139*IT_0188 + IT_0184*IT_0191 
      + IT_0026*IT_0192 + IT_0178*IT_0193 + IT_0167*IT_0194 + IT_0175*IT_0195);
    const ccomplex_t IT_0197 = cpow(s_12 + IT_0009 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0199 = 0.666666666666667*IT_0198;
    const ccomplex_t IT_0200 = cpow(IT_0199, 2);
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = IT_0197*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0204 = IT_0197*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = (-2)*IT_0203;
    const ccomplex_t IT_0207 = IT_0205*IT_0206;
    const ccomplex_t IT_0208 = 0.5*IT_0202 + (-0.0833333333333333)*IT_0207;
    const ccomplex_t IT_0209 = conj(IT_0139) + conj(IT_0167);
    const ccomplex_t IT_0210 = -s_14;
    const ccomplex_t IT_0211 = s_13 + IT_0210;
    const ccomplex_t IT_0212 = m_t*IT_0211;
    const ccomplex_t IT_0213 = 6*IT_0212;
    const ccomplex_t IT_0214 = 18*IT_0212;
    const ccomplex_t IT_0215 = (-18)*IT_0212;
    const ccomplex_t IT_0216 = (-6)*IT_0212;
    const ccomplex_t IT_0217 = s_34*IT_0009;
    const ccomplex_t IT_0218 = (-36)*IT_0217;
    const ccomplex_t IT_0219 = IT_0009*IT_0063;
    const ccomplex_t IT_0220 = (-36)*IT_0219;
    const ccomplex_t IT_0221 = s_13*s_14;
    const ccomplex_t IT_0222 = 72*IT_0221;
    const ccomplex_t IT_0223 = IT_0218 + IT_0220 + IT_0222;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = conj(IT_0062)*IT_0214;
    const ccomplex_t IT_0226 = (-18)*IT_0217;
    const ccomplex_t IT_0227 = (-18)*IT_0219;
    const ccomplex_t IT_0228 = 36*IT_0221;
    const ccomplex_t IT_0229 = IT_0226 + IT_0227 + IT_0228;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0208*IT_0230;
    const ccomplex_t IT_0232 = 0.25*IT_0207;
    const ccomplex_t IT_0233 = (-6)*IT_0217;
    const ccomplex_t IT_0234 = (-6)*IT_0219;
    const ccomplex_t IT_0235 = 12*IT_0221;
    const ccomplex_t IT_0236 = IT_0233 + IT_0234 + IT_0235;
    const ccomplex_t IT_0237 = 6*conj(IT_0023);
    const ccomplex_t IT_0238 = (-6)*IT_0023;
    const ccomplex_t IT_0239 = -IT_0236;
    const ccomplex_t IT_0240 = (-0.5)*conj(IT_0141)*IT_0213 + (-0.5)*IT_0209
      *IT_0214 + (-0.5)*IT_0168*IT_0215 + (-0.5)*IT_0141*IT_0216 + IT_0208*((-12
      )*IT_0217 + (-12)*IT_0219 + 24*IT_0221) + (-0.5)*(conj(IT_0175) + conj
      (IT_0184))*IT_0229 + (-0.5)*(IT_0175 + IT_0184)*IT_0230 + (-0.5)*IT_0224
      *IT_0232 + (-0.5)*(conj(IT_0026) + conj(IT_0178))*IT_0236 + (-0.5)*IT_0212
      *((-12)*IT_0062 + 12*conj(IT_0062) + IT_0237 + IT_0238) + (-0.5)*(IT_0026 
      + IT_0178)*IT_0239;
    const ccomplex_t IT_0241 = conj(IT_0062)*IT_0213;
    const ccomplex_t IT_0242 = IT_0208*IT_0239;
    const ccomplex_t IT_0243 = IT_0233 + IT_0235;
    const ccomplex_t IT_0244 = IT_0226 + IT_0228;
    const ccomplex_t IT_0245 = IT_0062*IT_0213;
    const ccomplex_t IT_0246 = IT_0208*IT_0236;
    const ccomplex_t IT_0247 = conj(IT_0175)*(IT_0141*IT_0176 + IT_0167
      *IT_0179 + IT_0139*IT_0182 + IT_0099*IT_0195 + IT_0184*IT_0227 + IT_0026
      *IT_0234 + IT_0181*IT_0238 + IT_0178*IT_0243 + IT_0175*IT_0244 + IT_0245 +
       IT_0246);
    const ccomplex_t IT_0248 = IT_0023*IT_0024 + 0.166666666666667*IT_0139
      *IT_0179 + 0.166666666666667*IT_0167*IT_0182 + 0.166666666666667*IT_0141
      *IT_0185 + 0.166666666666667*IT_0099*IT_0191 + 0.166666666666667*IT_0175
      *IT_0227 + 0.166666666666667*IT_0178*IT_0234 + 0.166666666666667*IT_0026
      *IT_0243 + 0.166666666666667*IT_0184*IT_0244 + 0.166666666666667*IT_0245 +
       0.166666666666667*IT_0246;
    const ccomplex_t IT_0249 = IT_0062*IT_0171;
    const ccomplex_t IT_0250 = IT_0208*IT_0229;
    const ccomplex_t IT_0251 = IT_0062*IT_0214;
    const ccomplex_t IT_0252 = 18*IT_0023*IT_0024 + IT_0139*IT_0176 + IT_0141
      *IT_0182 + IT_0167*IT_0185 + IT_0099*IT_0192 + IT_0178*IT_0227 + IT_0175
      *IT_0234 + IT_0184*IT_0243 + IT_0026*IT_0244 + IT_0250 + IT_0251;
    const ccomplex_t IT_0253 = conj(IT_0178)*(IT_0167*IT_0176 + IT_0141
      *IT_0179 + IT_0180*IT_0181 + IT_0139*IT_0185 + IT_0099*IT_0193 + IT_0026
      *IT_0227 + IT_0184*IT_0234 + IT_0175*IT_0243 + IT_0178*IT_0244 + IT_0250 +
       IT_0251);
    const ccomplex_t IT_0254 = conj(IT_0023)*(18*IT_0024*IT_0026 + IT_0023
      *IT_0065 + IT_0062*IT_0066 + IT_0099*IT_0101 + IT_0139*IT_0140) + conj
      (IT_0062)*(IT_0066*(IT_0023 + IT_0141) + IT_0062*((-36)*IT_0063 + IT_0142)
       + IT_0168*IT_0171) + IT_0063*(IT_0141*IT_0172 + IT_0167*IT_0173) +
       IT_0186 + IT_0196 + IT_0208*((18*conj(IT_0023) + (-36)*IT_0062 + 36*conj
      (IT_0062) + IT_0180)*IT_0212 + IT_0209*IT_0213 + conj(IT_0141)*IT_0214 +
       IT_0141*IT_0215 + IT_0168*IT_0216 + IT_0208*IT_0224) + IT_0026*(IT_0225 +
       IT_0231) + IT_0178*(IT_0172*IT_0181 + IT_0225 + IT_0231) + (-2)*IT_0232
      *IT_0240 + IT_0175*(IT_0173*IT_0181 + IT_0241 + IT_0242) + IT_0184*
      (IT_0024*IT_0237 + IT_0241 + IT_0242) + IT_0247 + 6*conj(IT_0184)*IT_0248 
      + conj(IT_0139)*(IT_0065*IT_0139 + IT_0023*IT_0140 + IT_0064*IT_0167 +
       IT_0026*IT_0176 + IT_0175*IT_0182 + IT_0141*IT_0183 + IT_0179*IT_0184 +
       IT_0178*IT_0185 + IT_0099*IT_0188 + IT_0249) + conj(IT_0167)*(IT_0064
      *IT_0139 + IT_0140*IT_0141 + IT_0065*IT_0167 + IT_0176*IT_0178 + IT_0175
      *IT_0179 + IT_0182*IT_0184 + IT_0026*IT_0185 + IT_0099*IT_0194 + IT_0063
      *IT_0238 + IT_0249) + conj(IT_0026)*IT_0252 + IT_0253;
    return create_ccomplex_return(IT_0254);
}

