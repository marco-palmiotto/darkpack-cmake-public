#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_anti_u_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_anti_u_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
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
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_00*conj(U_sb_01) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0005 = (-0.666666666666667)*IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0002;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0003;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0007 + (-3)*IT_0009);
    const ccomplex_t IT_0011 = (-0.166666666666667)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_sb_1, 2);
    const ccomplex_t IT_0013 = pow(m_sb_2, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0012 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0005*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = pow(m_u, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0013 + IT_0019 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = cpow(IT_0023, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *e_em*IT_0003*IT_0021*IT_0022*IT_0024*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0003*IT_0021*IT_0024*IT_0027*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0013 + IT_0019 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *e_em*IT_0003*IT_0021*IT_0022*IT_0024*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0003*IT_0021*IT_0024*IT_0027*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = IT_0033*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = -IT_0032 + -IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0043 = (-0.666666666666667)*IT_0042;
    const ccomplex_t IT_0044 = IT_0014*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0005*IT_0045;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = IT_0003*IT_0027;
    const ccomplex_t IT_0049 = e_em*conj(U_Wm2)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = cos(beta);
    const ccomplex_t IT_0052 = cpow(IT_0051, -1);
    const ccomplex_t IT_0053 = IT_0003*IT_0027*IT_0052;
    const ccomplex_t IT_0054 = m_b*conj(U_d2)*e_em*IT_0021*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0050 + (-0.5)*IT_0056);
    const ccomplex_t IT_0058 = IT_0003*IT_0022;
    const ccomplex_t IT_0059 = e_em*U_Wm2*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = IT_0003*IT_0022*IT_0052;
    const ccomplex_t IT_0062 = m_b*U_d2*e_em*IT_0021*conj(U_sb_11)*V_ub_mod;
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0060 + (-0.5)*IT_0064);
    const ccomplex_t IT_0066 = IT_0057*IT_0065;
    const ccomplex_t IT_0067 = IT_0033*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = e_em*conj(U_Wm1)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0070 = IT_0048*IT_0069;
    const ccomplex_t IT_0071 = m_b*conj(U_d1)*e_em*IT_0021*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0072 = IT_0053*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0070 + (-0.5)*IT_0073);
    const ccomplex_t IT_0075 = e_em*U_Wm1*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0076 = IT_0058*IT_0075;
    const ccomplex_t IT_0077 = m_b*U_d1*e_em*IT_0021*conj(U_sb_11)*V_ub_mod;
    const ccomplex_t IT_0078 = IT_0061*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0076 + (-0.5)*IT_0079);
    const ccomplex_t IT_0081 = IT_0074*IT_0080;
    const ccomplex_t IT_0082 = IT_0020*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = -IT_0068 + -IT_0083;
    const ccomplex_t IT_0085 = IT_0012*IT_0019;
    const ccomplex_t IT_0086 = (-0.333333333333333)*IT_0004;
    const ccomplex_t IT_0087 = m_u*IT_0086;
    const ccomplex_t IT_0088 = m_u*IT_0005;
    const ccomplex_t IT_0089 = IT_0016*(IT_0087 + -IT_0088);
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = IT_0045*(IT_0087 + -IT_0088);
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = sin(alpha);
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*m_u*e_em*IT_0003*IT_0021
      *IT_0024*IT_0093;
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = cpow(V_tb, 2);
    const ccomplex_t IT_0097 = cos(alpha);
    const ccomplex_t IT_0098 = m_b*mu_h*IT_0021;
    const ccomplex_t IT_0099 = -IT_0023*IT_0093;
    const ccomplex_t IT_0100 = IT_0051*IT_0097;
    const ccomplex_t IT_0101 = IT_0099 + IT_0100;
    const ccomplex_t IT_0102 = cpow(IT_0000, -2);
    const ccomplex_t IT_0103 = pow(m_b, 2);
    const ccomplex_t IT_0104 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0101*(IT_0003 + 1.
      /3*IT_0002*IT_0102) + IT_0003*IT_0021*IT_0052*IT_0097*IT_0103)*conj
      (U_sb_01) + -1./2*IT_0003*IT_0052*IT_0093*IT_0098*conj(U_sb_11)) + U_sb_10
      *(IT_0003*IT_0052*IT_0093*IT_0098*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0002*IT_0101*IT_0102 + (-3)*IT_0003*IT_0021*IT_0052*IT_0097*IT_0103)
      *conj(U_sb_11));
    const ccomplex_t IT_0105 = (0 + _Complex_I*-1)*(e_em*IT_0104 + (
      -1.4142135623731)*A_b*IT_0096*IT_0097*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = IT_0095*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*m_u*e_em*IT_0003*IT_0021
      *IT_0024*IT_0097;
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = IT_0051*IT_0093;
    const ccomplex_t IT_0114 = IT_0023*IT_0097;
    const ccomplex_t IT_0115 = IT_0113 + IT_0114;
    const ccomplex_t IT_0116 = 2*U_sb_00*((IT_0003*IT_0021*IT_0052*IT_0093
      *IT_0103 + (-0.5)*m_W*(IT_0003 + 1./3*IT_0002*IT_0102)*IT_0115)*conj
      (U_sb_01) + 1./2*IT_0003*IT_0052*IT_0097*IT_0098*conj(U_sb_11)) + U_sb_10*
      (IT_0003*IT_0052*IT_0097*IT_0098*conj(U_sb_01) + 2*(IT_0003*IT_0021
      *IT_0052*IT_0093*IT_0103 + -1./3*m_W*IT_0002*IT_0102*IT_0115)*conj(U_sb_11
      ));
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*(e_em*IT_0116 +
       1.4142135623731*A_b*IT_0093*IT_0096*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = IT_0112*IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0012 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0121 = IT_0119*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = -IT_0110 + -IT_0122;
    const ccomplex_t IT_0124 = 18*s_34;
    const ccomplex_t IT_0125 = (-18)*IT_0019;
    const ccomplex_t IT_0126 = IT_0124 + IT_0125;
    const ccomplex_t IT_0127 = s_34*IT_0012;
    const ccomplex_t IT_0128 = (-6)*IT_0127;
    const ccomplex_t IT_0129 = s_13*s_14;
    const ccomplex_t IT_0130 = 12*IT_0129;
    const ccomplex_t IT_0131 = IT_0128 + IT_0130;
    const ccomplex_t IT_0132 = (-18)*IT_0127;
    const ccomplex_t IT_0133 = 36*IT_0129;
    const ccomplex_t IT_0134 = IT_0132 + IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0136 = IT_0029*IT_0080;
    const ccomplex_t IT_0137 = IT_0020*IT_0135*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0139 = IT_0037*IT_0065;
    const ccomplex_t IT_0140 = IT_0033*IT_0138*IT_0139;
    const ccomplex_t IT_0141 = -IT_0137 + -IT_0140;
    const ccomplex_t IT_0142 = 6*s_34;
    const ccomplex_t IT_0143 = IT_0090*IT_0142;
    const ccomplex_t IT_0144 = IT_0124*IT_0141;
    const ccomplex_t IT_0145 = IT_0143 + IT_0144;
    const ccomplex_t IT_0146 = IT_0026*IT_0074;
    const ccomplex_t IT_0147 = IT_0020*IT_0135*IT_0146;
    const ccomplex_t IT_0148 = IT_0035*IT_0057;
    const ccomplex_t IT_0149 = IT_0033*IT_0138*IT_0148;
    const ccomplex_t IT_0150 = m_u*IT_0032;
    const ccomplex_t IT_0151 = m_u*IT_0040;
    const ccomplex_t IT_0152 = -IT_0147 + -IT_0149 + IT_0150 + IT_0151;
    const ccomplex_t IT_0153 = IT_0089 + IT_0092;
    const ccomplex_t IT_0154 = conj(IT_0089) + conj(IT_0092);
    const ccomplex_t IT_0155 = conj(IT_0090) + conj(IT_0091);
    const ccomplex_t IT_0156 = IT_0068 + IT_0083;
    const ccomplex_t IT_0157 = m_u*IT_0156;
    const ccomplex_t IT_0158 = conj(IT_0141) + conj(IT_0157);
    const ccomplex_t IT_0159 = -IT_0019;
    const ccomplex_t IT_0160 = 6*IT_0123;
    const ccomplex_t IT_0161 = 6*conj(IT_0123);
    const ccomplex_t IT_0162 = IT_0090*IT_0124;
    const ccomplex_t IT_0163 = -s_14;
    const ccomplex_t IT_0164 = s_13 + IT_0163;
    const ccomplex_t IT_0165 = m_u*IT_0164;
    const ccomplex_t IT_0166 = 18*IT_0123;
    const ccomplex_t IT_0167 = 18*conj(IT_0123);
    const ccomplex_t IT_0168 = s_34 + IT_0019;
    const ccomplex_t IT_0169 = m_u*e_em*IT_0003*IT_0021;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = m_b*e_em*mu_h*IT_0021*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0172 = IT_0003*IT_0023*IT_0052;
    const ccomplex_t IT_0173 = IT_0171*IT_0172;
    const ccomplex_t IT_0174 = m_b*e_em*mu_h*IT_0021*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0175 = IT_0172*IT_0174;
    const ccomplex_t IT_0176 = A_b*IT_0096*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0177 = IT_0051*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = A_b*IT_0096*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0180 = IT_0051*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = (-0.5)*IT_0173 + 0.5*IT_0175 + 0.5*IT_0178 + (
      -0.5)*IT_0181;
    const ccomplex_t IT_0183 = IT_0170*IT_0182;
    const ccomplex_t IT_0184 = pow(m_Z, 2);
    const ccomplex_t IT_0185 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + IT_0184 
      + -reg_prop, -1);
    const ccomplex_t IT_0186 = IT_0183*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0189 = IT_0003*IT_0171;
    const ccomplex_t IT_0190 = IT_0003*IT_0174;
    const ccomplex_t IT_0191 = IT_0023*IT_0176;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = IT_0023*IT_0179;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = (-0.5)*IT_0189 + 0.5*IT_0190 + (-0.5)*IT_0192 +
       0.5*IT_0194;
    const ccomplex_t IT_0196 = m_u*e_em*IT_0003*IT_0021*IT_0024*IT_0051;
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = IT_0195*IT_0197;
    const ccomplex_t IT_0199 = IT_0188*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = -IT_0187 + -IT_0200;
    const ccomplex_t IT_0202 = (-6)*IT_0201;
    const ccomplex_t IT_0203 = (-6)*conj(IT_0201);
    const ccomplex_t IT_0204 = s_13 + s_14;
    const ccomplex_t IT_0205 = m_u*IT_0204;
    const ccomplex_t IT_0206 = (-18)*IT_0152;
    const ccomplex_t IT_0207 = (-6)*IT_0089;
    const ccomplex_t IT_0208 = (-6)*IT_0092;
    const ccomplex_t IT_0209 = IT_0206 + IT_0207 + IT_0208;
    const ccomplex_t IT_0210 = (-18)*conj(IT_0152);
    const ccomplex_t IT_0211 = (-6)*conj(IT_0089);
    const ccomplex_t IT_0212 = (-6)*conj(IT_0092);
    const ccomplex_t IT_0213 = IT_0210 + IT_0211 + IT_0212;
    const ccomplex_t IT_0214 = (-18)*IT_0089;
    const ccomplex_t IT_0215 = IT_0090 + IT_0091;
    const ccomplex_t IT_0216 = (-18)*conj(IT_0089);
    const ccomplex_t IT_0217 = (-18)*IT_0090;
    const ccomplex_t IT_0218 = (-18)*IT_0091;
    const ccomplex_t IT_0219 = (-6)*IT_0141;
    const ccomplex_t IT_0220 = (-18)*conj(IT_0090);
    const ccomplex_t IT_0221 = (-18)*conj(IT_0091);
    const ccomplex_t IT_0222 = (-6)*conj(IT_0141);
    const ccomplex_t IT_0223 = (-18)*IT_0141;
    const ccomplex_t IT_0224 = (-6)*IT_0090;
    const ccomplex_t IT_0225 = (-6)*IT_0091;
    const ccomplex_t IT_0226 = (-18)*conj(IT_0141);
    const ccomplex_t IT_0227 = (-6)*conj(IT_0090);
    const ccomplex_t IT_0228 = (-6)*conj(IT_0091);
    const ccomplex_t IT_0229 = conj(IT_0141)*IT_0207 + conj(IT_0157)*IT_0209 +
       IT_0141*IT_0211 + IT_0157*IT_0213 + IT_0155*IT_0214 + IT_0215*IT_0216 +
       conj(IT_0092)*(IT_0217 + IT_0218 + IT_0219) + IT_0092*(IT_0220 + IT_0221 
      + IT_0222) + conj(IT_0152)*(IT_0223 + IT_0224 + IT_0225) + IT_0152*
      (IT_0226 + IT_0227 + IT_0228);
    const ccomplex_t IT_0230 = m_u*s_14;
    const ccomplex_t IT_0231 = conj(IT_0084)*IT_0209 + IT_0084*IT_0213 + (-18)
      *conj(IT_0018)*(IT_0092 + 1./3*IT_0152 + (-0.0555555555555555)*IT_0214) + 
      (-18)*IT_0018*(conj(IT_0092) + 1./3*conj(IT_0152) + (-0.0555555555555555)
      *IT_0216) + (-6)*conj(IT_0047)*(IT_0157 + (-0.166666666666667)*IT_0217 + (
      -0.166666666666667)*IT_0218 + (-0.166666666666667)*IT_0219) + (-6)*IT_0047
      *(conj(IT_0157) + (-0.166666666666667)*IT_0220 + (-0.166666666666667)
      *IT_0221 + (-0.166666666666667)*IT_0222) + (-18)*conj(IT_0041)*(IT_0157 + 
      (-0.0555555555555556)*IT_0223 + (-0.0555555555555556)*IT_0224 + (
      -0.0555555555555556)*IT_0225) + (-18)*IT_0041*(conj(IT_0157) + (
      -0.0555555555555556)*IT_0226 + (-0.0555555555555556)*IT_0227 + (
      -0.0555555555555556)*IT_0228);
    const ccomplex_t IT_0232 = 36*s_34;
    const ccomplex_t IT_0233 = 18*IT_0019;
    const ccomplex_t IT_0234 = IT_0124 + IT_0233;
    const ccomplex_t IT_0235 = (-18)*s_34;
    const ccomplex_t IT_0236 = IT_0125 + IT_0235;
    const ccomplex_t IT_0237 = 6*IT_0152;
    const ccomplex_t IT_0238 = 6*IT_0041;
    const ccomplex_t IT_0239 = 6*conj(IT_0152);
    const ccomplex_t IT_0240 = 6*conj(IT_0041);
    const ccomplex_t IT_0241 = m_u*s_13;
    const ccomplex_t IT_0242 = conj(IT_0084)*(IT_0090 + IT_0091 + 3*IT_0141 +
       3*IT_0157) + 0.166666666666667*conj(IT_0018)*(18*IT_0090 + 18*IT_0091 + 6
      *IT_0141 + 6*IT_0157) + IT_0084*(conj(IT_0090) + conj(IT_0091) + 3*conj
      (IT_0141) + 3*conj(IT_0157)) + 0.166666666666667*IT_0018*(18*conj(IT_0090)
       + 18*conj(IT_0091) + 6*conj(IT_0141) + 6*conj(IT_0157)) +
       0.166666666666667*conj(IT_0047)*(18*IT_0089 + 18*IT_0092 + IT_0237) +
       0.166666666666667*(conj(IT_0089) + conj(IT_0092) + 3*conj(IT_0152))
      *IT_0238 + 0.166666666666667*IT_0047*(18*conj(IT_0089) + 18*conj(IT_0092) 
      + IT_0239) + 0.166666666666667*(IT_0089 + IT_0092 + 3*IT_0152)*IT_0240;
    const ccomplex_t IT_0243 = (-6)*(conj(IT_0018)*(IT_0041 + 3*IT_0047) +
       IT_0018*(conj(IT_0041) + 3*conj(IT_0047)) + 3*(conj(IT_0041) +
       0.333333333333333*conj(IT_0047))*IT_0084 + 3*(IT_0041 + 0.333333333333333
      *IT_0047)*conj(IT_0084))*IT_0085 + ((conj(IT_0089) + conj(IT_0090) + conj
      (IT_0091) + conj(IT_0092))*IT_0123 + (IT_0089 + IT_0090 + IT_0091 +
       IT_0092)*conj(IT_0123))*IT_0126 + (conj(IT_0041)*IT_0047 + IT_0041*conj
      (IT_0047) + conj(IT_0018)*IT_0084 + IT_0018*conj(IT_0084))*IT_0131 + 
      (IT_0018*conj(IT_0018) + IT_0041*conj(IT_0041) + IT_0047*conj(IT_0047) +
       IT_0084*conj(IT_0084))*IT_0134 + conj(IT_0141)*IT_0145 + IT_0142*(conj
      (IT_0152)*IT_0153 + IT_0152*IT_0154 + IT_0155*(IT_0141 + IT_0157) +
       IT_0091*IT_0158) + IT_0124*(IT_0152*conj(IT_0152) + IT_0153*IT_0154 +
       IT_0091*IT_0155 + IT_0157*IT_0158) + (s_34 + IT_0159)*((conj(IT_0141) +
       conj(IT_0152) + conj(IT_0157))*IT_0160 + (IT_0141 + IT_0152 + IT_0157)
      *IT_0161) + IT_0155*IT_0162 + IT_0165*((conj(IT_0041) + conj(IT_0084))
      *IT_0160 + (IT_0041 + IT_0084)*IT_0161 + (conj(IT_0018) + conj(IT_0047))
      *IT_0166 + (IT_0018 + IT_0047)*IT_0167) + conj(IT_0157)*(IT_0145 + IT_0168
      *IT_0202) + IT_0157*IT_0168*IT_0203 + (conj(IT_0084)*IT_0202 + IT_0084
      *IT_0203)*IT_0205 + IT_0019*IT_0229 + IT_0230*IT_0231 + IT_0123*conj
      (IT_0123)*((-36)*IT_0019 + IT_0232) + conj(IT_0201)*(IT_0153*IT_0234 +
       IT_0215*IT_0236 + IT_0168*(IT_0219 + IT_0237) + IT_0205*((-18)*IT_0018 +
       18*IT_0047 + IT_0238)) + IT_0201*(conj(IT_0201)*(36*IT_0019 + IT_0232) +
       IT_0154*IT_0234 + IT_0155*IT_0236 + IT_0168*(IT_0222 + IT_0239) + IT_0205
      *((-18)*conj(IT_0018) + 18*conj(IT_0047) + IT_0240)) + 6*IT_0241*IT_0242;
    return create_ccomplex_return(IT_0243);
}

