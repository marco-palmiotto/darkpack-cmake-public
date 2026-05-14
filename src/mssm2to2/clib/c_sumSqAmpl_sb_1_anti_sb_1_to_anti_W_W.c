#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_anti_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t V_cb = param->V_cb;
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t IT_0000 = pow(m_sb_1, 2);
    const ccomplex_t IT_0001 = cpow(s_12 + IT_0000 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*IT_0001;
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0004 = -IT_0003;
    const ccomplex_t IT_0005 = 0.333333333333333*IT_0003;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 2*IT_0003;
    const ccomplex_t IT_0008 = IT_0005*IT_0007;
    const ccomplex_t IT_0009 = IT_0006 + IT_0008;
    const ccomplex_t IT_0010 = IT_0002*IT_0009;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0012*conj(U_sb_00);
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0012*U_sb_00;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0000 + IT_0018 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0000 + IT_0018 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0012
      *IT_0024*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0012
      *IT_0027*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = IT_0023*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_13 + IT_0000 + IT_0018 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0012*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0012*U_sb_00*conj(U_st_00);
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = IT_0033*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0000 + IT_0018 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0012*conj(U_sb_00)*U_st_01;
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0012*U_sb_00*conj(U_st_01);
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = IT_0041*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = cos(theta_W);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*e_em*IT_0012*IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = cpow(IT_0050, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0012
      *IT_0050 + 1./3*IT_0011*IT_0053)*U_sb_00*conj(U_sb_00) + IT_0011*IT_0053
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0055 = (-0.333333333333333)*IT_0054;
    const ccomplex_t IT_0056 = IT_0052*IT_0055;
    const ccomplex_t IT_0057 = 2*IT_0051;
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = IT_0049*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0010 + IT_0021 + IT_0031 + IT_0039 +
       IT_0047 + -IT_0060;
    const ccomplex_t IT_0062 = pow(m_W, -2);
    const ccomplex_t IT_0063 = pow(s_13, 2);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = -IT_0000;
    const ccomplex_t IT_0066 = pow(s_14, 2);
    const ccomplex_t IT_0067 = IT_0062*IT_0066;
    const ccomplex_t IT_0068 = (IT_0064 + IT_0065)*(IT_0065 + IT_0067);
    const ccomplex_t IT_0069 = 0.666666666666667*IT_0003;
    const ccomplex_t IT_0070 = IT_0007*IT_0069;
    const ccomplex_t IT_0071 = IT_0004*IT_0069;
    const ccomplex_t IT_0072 = -IT_0070 + -IT_0071;
    const ccomplex_t IT_0073 = IT_0009 + IT_0072;
    const ccomplex_t IT_0074 = IT_0002*IT_0073;
    const ccomplex_t IT_0075 = (-0.666666666666667)*IT_0054;
    const ccomplex_t IT_0076 = IT_0052*IT_0075;
    const ccomplex_t IT_0077 = IT_0057*IT_0075;
    const ccomplex_t IT_0078 = -IT_0076 + -IT_0077;
    const ccomplex_t IT_0079 = IT_0059 + IT_0078;
    const ccomplex_t IT_0080 = IT_0049*IT_0079;
    const ccomplex_t IT_0081 = 0.5*IT_0074 + -IT_0080;
    const ccomplex_t IT_0082 = s_13*s_23*IT_0062;
    const ccomplex_t IT_0083 = s_14*s_24*IT_0062;
    const ccomplex_t IT_0084 = (s_12 + -IT_0082)*(s_12 + -IT_0083);
    const ccomplex_t IT_0085 = 3*IT_0084;
    const ccomplex_t IT_0086 = IT_0003*IT_0005;
    const ccomplex_t IT_0087 = IT_0003*IT_0069;
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = IT_0086 + IT_0088;
    const ccomplex_t IT_0090 = IT_0002*IT_0089;
    const ccomplex_t IT_0091 = IT_0051*IT_0055;
    const ccomplex_t IT_0092 = IT_0051*IT_0075;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = IT_0091 + IT_0093;
    const ccomplex_t IT_0095 = IT_0049*IT_0094;
    const ccomplex_t IT_0096 = (-0.5)*IT_0090 + IT_0095;
    const ccomplex_t IT_0097 = s_14*s_34*IT_0062;
    const ccomplex_t IT_0098 = (s_12 + -IT_0082)*(s_13 + -IT_0097);
    const ccomplex_t IT_0099 = 3*IT_0098;
    const ccomplex_t IT_0100 = s_12 + IT_0000;
    const ccomplex_t IT_0101 = IT_0071*IT_0100;
    const ccomplex_t IT_0102 = (-2)*IT_0003;
    const ccomplex_t IT_0103 = IT_0005*IT_0102;
    const ccomplex_t IT_0104 = s_13 + s_23;
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = IT_0101 + IT_0105;
    const ccomplex_t IT_0107 = IT_0069*IT_0102;
    const ccomplex_t IT_0108 = s_23*IT_0107;
    const ccomplex_t IT_0109 = IT_0006*IT_0100;
    const ccomplex_t IT_0110 = -IT_0108 + (-2)*IT_0109;
    const ccomplex_t IT_0111 = IT_0106 + IT_0110;
    const ccomplex_t IT_0112 = IT_0002*IT_0111;
    const ccomplex_t IT_0113 = IT_0076*IT_0100;
    const ccomplex_t IT_0114 = (-2)*IT_0051;
    const ccomplex_t IT_0115 = IT_0055*IT_0114;
    const ccomplex_t IT_0116 = IT_0104*IT_0115;
    const ccomplex_t IT_0117 = IT_0113 + IT_0116;
    const ccomplex_t IT_0118 = IT_0075*IT_0114;
    const ccomplex_t IT_0119 = s_23*IT_0118;
    const ccomplex_t IT_0120 = IT_0056*IT_0100;
    const ccomplex_t IT_0121 = -IT_0119 + (-2)*IT_0120;
    const ccomplex_t IT_0122 = IT_0117 + IT_0121;
    const ccomplex_t IT_0123 = IT_0049*IT_0122;
    const ccomplex_t IT_0124 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0125 = cos(alpha);
    const ccomplex_t IT_0126 = sin(beta);
    const ccomplex_t IT_0127 = cos(beta);
    const ccomplex_t IT_0128 = sin(alpha);
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_W*e_em*IT_0012*(IT_0125
      *IT_0126 + -IT_0127*IT_0128);
    const ccomplex_t IT_0130 = cpow(V_tb, 2);
    const ccomplex_t IT_0131 = pow(m_W, -1);
    const ccomplex_t IT_0132 = m_b*mu_h*IT_0131;
    const ccomplex_t IT_0133 = cpow(IT_0127, -1);
    const ccomplex_t IT_0134 = IT_0127*IT_0128;
    const ccomplex_t IT_0135 = IT_0125*IT_0126;
    const ccomplex_t IT_0136 = IT_0134 + IT_0135;
    const ccomplex_t IT_0137 = cpow(IT_0050, -2);
    const ccomplex_t IT_0138 = pow(m_b, 2);
    const ccomplex_t IT_0139 = 2*U_sb_00*(((-0.5)*m_W*IT_0136*(IT_0012 + 1./3
      *IT_0011*IT_0137) + IT_0012*IT_0128*IT_0131*IT_0133*IT_0138)*conj(U_sb_00)
       + 1./2*IT_0012*IT_0125*IT_0132*IT_0133*conj(U_sb_10)) + U_sb_10*(IT_0012
      *IT_0125*IT_0132*IT_0133*conj(U_sb_00) + (-0.666666666666667)*(m_W*IT_0011
      *IT_0136*IT_0137 + (-3)*IT_0012*IT_0128*IT_0131*IT_0133*IT_0138)*conj
      (U_sb_10));
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(e_em*IT_0139 +
       1.4142135623731*A_b*IT_0128*IT_0130*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0129*IT_0141;
    const ccomplex_t IT_0143 = IT_0124*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_W*e_em*IT_0012*(IT_0125
      *IT_0127 + IT_0126*IT_0128);
    const ccomplex_t IT_0147 = -IT_0126*IT_0128;
    const ccomplex_t IT_0148 = IT_0125*IT_0127;
    const ccomplex_t IT_0149 = IT_0147 + IT_0148;
    const ccomplex_t IT_0150 = (-2)*U_sb_00*((IT_0012*IT_0125*IT_0131*IT_0133
      *IT_0138 + (-0.5)*m_W*(IT_0012 + 1./3*IT_0011*IT_0137)*IT_0149)*conj
      (U_sb_00) + -1./2*IT_0012*IT_0128*IT_0132*IT_0133*conj(U_sb_10)) + U_sb_10
      *(IT_0012*IT_0128*IT_0132*IT_0133*conj(U_sb_00) + (-2)*(IT_0012*IT_0125
      *IT_0131*IT_0133*IT_0138 + -1./3*m_W*IT_0011*IT_0137*IT_0149)*conj(U_sb_10
      ));
    const ccomplex_t IT_0151 = (0 + _Complex_I*-1)*(e_em*IT_0150 + (
      -1.4142135623731)*A_b*IT_0125*IT_0130*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = IT_0146*IT_0152;
    const ccomplex_t IT_0154 = IT_0145*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = pow(e_em, 2);
    const ccomplex_t IT_0157 = cpow(IT_0011, -2);
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0156*IT_0157*U_sb_00*conj
      (U_sb_00);
    const ccomplex_t IT_0159 = 0.5*IT_0158;
    const ccomplex_t IT_0160 = (-0.5)*IT_0112 + IT_0123 + -IT_0144 + -IT_0155 
      + IT_0159;
    const ccomplex_t IT_0161 = pow(m_W, -4);
    const ccomplex_t IT_0162 = s_13*s_14*s_34*IT_0161;
    const ccomplex_t IT_0163 = IT_0064 + IT_0065 + IT_0067 + -IT_0162;
    const ccomplex_t IT_0164 = (-3)*IT_0163;
    const ccomplex_t IT_0165 = -IT_0034;
    const ccomplex_t IT_0166 = IT_0037*IT_0165;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0038 + IT_0167;
    const ccomplex_t IT_0169 = IT_0033*IT_0168;
    const ccomplex_t IT_0170 = -IT_0070;
    const ccomplex_t IT_0171 = IT_0009 + IT_0170;
    const ccomplex_t IT_0172 = IT_0002*IT_0171;
    const ccomplex_t IT_0173 = -IT_0013;
    const ccomplex_t IT_0174 = IT_0016*IT_0173;
    const ccomplex_t IT_0175 = -IT_0174;
    const ccomplex_t IT_0176 = IT_0017 + IT_0175;
    const ccomplex_t IT_0177 = IT_0020*IT_0176;
    const ccomplex_t IT_0178 = -IT_0025;
    const ccomplex_t IT_0179 = IT_0029*IT_0178;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = IT_0030 + IT_0180;
    const ccomplex_t IT_0182 = IT_0023*IT_0181;
    const ccomplex_t IT_0183 = -IT_0042;
    const ccomplex_t IT_0184 = IT_0045*IT_0183;
    const ccomplex_t IT_0185 = -IT_0184;
    const ccomplex_t IT_0186 = IT_0046 + IT_0185;
    const ccomplex_t IT_0187 = IT_0041*IT_0186;
    const ccomplex_t IT_0188 = -IT_0077;
    const ccomplex_t IT_0189 = IT_0059 + IT_0188;
    const ccomplex_t IT_0190 = IT_0049*IT_0189;
    const ccomplex_t IT_0191 = IT_0169 + 0.5*IT_0172 + IT_0177 + IT_0182 +
       IT_0187 + -IT_0190;
    const ccomplex_t IT_0192 = -(IT_0000 + -IT_0064)*(s_12 + -IT_0083);
    const ccomplex_t IT_0193 = (-3)*IT_0192;
    const ccomplex_t IT_0194 = IT_0049*IT_0091;
    const ccomplex_t IT_0195 = IT_0002*IT_0086;
    const ccomplex_t IT_0196 = -IT_0021 + -IT_0031 + -IT_0039 + -IT_0047 +
       IT_0194 + (-0.5)*IT_0195;
    const ccomplex_t IT_0197 = -(IT_0000 + -IT_0064)*(s_13 + -IT_0097);
    const ccomplex_t IT_0198 = (-3)*IT_0197;
    const ccomplex_t IT_0199 = -IT_0071;
    const ccomplex_t IT_0200 = IT_0009 + IT_0199;
    const ccomplex_t IT_0201 = IT_0002*IT_0200;
    const ccomplex_t IT_0202 = -IT_0076;
    const ccomplex_t IT_0203 = IT_0059 + IT_0202;
    const ccomplex_t IT_0204 = IT_0049*IT_0203;
    const ccomplex_t IT_0205 = 0.5*IT_0201 + -IT_0204;
    const ccomplex_t IT_0206 = -(IT_0000 + -IT_0067)*(s_12 + -IT_0082);
    const ccomplex_t IT_0207 = (-3)*IT_0206;
    const ccomplex_t IT_0208 = 3*IT_0061*(conj(IT_0061)*IT_0068 +
       0.333333333333333*conj(IT_0081)*IT_0085 + 0.333333333333333*conj(IT_0096)
      *IT_0099 + 0.333333333333333*conj(IT_0160)*IT_0164 + 0.333333333333333
      *conj(IT_0191)*IT_0193 + 0.333333333333333*conj(IT_0196)*IT_0198 +
       0.333333333333333*conj(IT_0205)*IT_0207);
    const ccomplex_t IT_0209 = pow(s_24, 2);
    const ccomplex_t IT_0210 = IT_0062*IT_0209;
    const ccomplex_t IT_0211 = -(IT_0000 + -IT_0064)*(IT_0000 + -IT_0210);
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = s_13*s_24*s_34*IT_0161;
    const ccomplex_t IT_0214 = -IT_0083;
    const ccomplex_t IT_0215 = -IT_0082;
    const ccomplex_t IT_0216 = s_12 + IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = IT_0213 + IT_0216;
    const ccomplex_t IT_0218 = 3*IT_0217;
    const ccomplex_t IT_0219 = s_24*s_34*IT_0062;
    const ccomplex_t IT_0220 = (s_12 + -IT_0082)*(s_23 + -IT_0219);
    const ccomplex_t IT_0221 = 3*IT_0220;
    const ccomplex_t IT_0222 = -(IT_0000 + -IT_0064)*(s_23 + -IT_0219);
    const ccomplex_t IT_0223 = (-3)*IT_0222;
    const ccomplex_t IT_0224 = (s_12 + -IT_0082)*(IT_0000 + -IT_0210);
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = (-3)*IT_0225;
    const ccomplex_t IT_0227 = 3*IT_0191*(0.333333333333333*conj(IT_0061)
      *IT_0193 + 0.333333333333333*IT_0085*conj(IT_0205) + conj(IT_0191)*IT_0212
       + 0.333333333333333*conj(IT_0160)*IT_0218 + 0.333333333333333*conj
      (IT_0096)*IT_0221 + 0.333333333333333*conj(IT_0196)*IT_0223 +
       0.333333333333333*conj(IT_0081)*IT_0226);
    const ccomplex_t IT_0228 = pow(s_23, 2);
    const ccomplex_t IT_0229 = IT_0062*IT_0228;
    const ccomplex_t IT_0230 = -(IT_0000 + -IT_0067)*(IT_0000 + -IT_0229);
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = s_14*s_23*s_34*IT_0161;
    const ccomplex_t IT_0233 = IT_0216 + IT_0232;
    const ccomplex_t IT_0234 = 3*IT_0233;
    const ccomplex_t IT_0235 = (s_12 + -IT_0083)*(IT_0000 + -IT_0229);
    const ccomplex_t IT_0236 = -IT_0235;
    const ccomplex_t IT_0237 = (-3)*IT_0236;
    const ccomplex_t IT_0238 = (s_13 + -IT_0097)*(IT_0000 + -IT_0229);
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = (-3)*IT_0239;
    const ccomplex_t IT_0241 = 3*IT_0205*(0.333333333333333*IT_0085*conj
      (IT_0191) + 0.333333333333333*IT_0099*conj(IT_0196) + 0.333333333333333
      *conj(IT_0061)*IT_0207 + conj(IT_0205)*IT_0231 + 0.333333333333333*conj
      (IT_0160)*IT_0234 + 0.333333333333333*conj(IT_0081)*IT_0237 +
       0.333333333333333*conj(IT_0096)*IT_0240);
    const ccomplex_t IT_0242 = (IT_0000 + -IT_0210)*(IT_0000 + -IT_0229);
    const ccomplex_t IT_0243 = s_23*s_24*s_34*IT_0161;
    const ccomplex_t IT_0244 = -IT_0210;
    const ccomplex_t IT_0245 = -IT_0229;
    const ccomplex_t IT_0246 = IT_0000 + IT_0243 + IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = 3*IT_0246;
    const ccomplex_t IT_0248 = (s_23 + -IT_0219)*(IT_0000 + -IT_0229);
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = (-3)*IT_0249;
    const ccomplex_t IT_0251 = 3*IT_0081*(0.333333333333333*conj(IT_0061)
      *IT_0085 + 0.333333333333333*conj(IT_0196)*IT_0221 + 0.333333333333333
      *conj(IT_0191)*IT_0226 + 0.333333333333333*conj(IT_0205)*IT_0237 + conj
      (IT_0081)*IT_0242 + 0.333333333333333*conj(IT_0160)*IT_0247 +
       0.333333333333333*conj(IT_0096)*IT_0250);
    const ccomplex_t IT_0252 = pow(s_34, 2);
    const ccomplex_t IT_0253 = IT_0161*IT_0252;
    const ccomplex_t IT_0254 = s_34*IT_0062;
    const ccomplex_t IT_0255 = s_13*s_34*IT_0062;
    const ccomplex_t IT_0256 = -IT_0255;
    const ccomplex_t IT_0257 = s_14 + IT_0256;
    const ccomplex_t IT_0258 = IT_0254*IT_0257;
    const ccomplex_t IT_0259 = (-3)*IT_0258;
    const ccomplex_t IT_0260 = s_23*s_34*IT_0062;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = s_24 + IT_0261;
    const ccomplex_t IT_0263 = IT_0254*IT_0262;
    const ccomplex_t IT_0264 = (-3)*IT_0263;
    const ccomplex_t IT_0265 = 0.166666666666667*conj(IT_0061)*IT_0164 +
       0.166666666666667*conj(IT_0191)*IT_0218 + 0.166666666666667*conj(IT_0205)
      *IT_0234 + 0.166666666666667*conj(IT_0081)*IT_0247 + conj(IT_0160)*(1 +
       0.5*IT_0253) + 0.166666666666667*conj(IT_0196)*IT_0259 +
       0.166666666666667*conj(IT_0096)*IT_0264;
    const ccomplex_t IT_0266 = IT_0062*IT_0252;
    const ccomplex_t IT_0267 = -(IT_0000 + -IT_0064)*(IT_0018 + -IT_0266);
    const ccomplex_t IT_0268 = (s_12 + -IT_0082)*(IT_0018 + -IT_0266);
    const ccomplex_t IT_0269 = 3*IT_0268;
    const ccomplex_t IT_0270 = (-0.333333333333333)*conj(IT_0061);
    const ccomplex_t IT_0271 = (-0.333333333333333)*conj(IT_0191);
    const ccomplex_t IT_0272 = (-0.333333333333333)*conj(IT_0205);
    const ccomplex_t IT_0273 = (-0.333333333333333)*conj(IT_0081);
    const ccomplex_t IT_0274 = (IT_0000 + -IT_0229)*(IT_0018 + -IT_0266);
    const ccomplex_t IT_0275 = -IT_0274;
    const ccomplex_t IT_0276 = IT_0208 + IT_0227 + IT_0241 + IT_0251 + 6
      *IT_0160*IT_0265 + (-3)*IT_0196*((-0.333333333333333)*conj(IT_0160)
      *IT_0259 + conj(IT_0196)*IT_0267 + (-0.333333333333333)*conj(IT_0096)
      *IT_0269 + IT_0198*IT_0270 + IT_0223*IT_0271 + IT_0099*IT_0272 + IT_0221
      *IT_0273) + (-3)*IT_0096*((-0.333333333333333)*conj(IT_0160)*IT_0264 + (
      -0.333333333333333)*conj(IT_0196)*IT_0269 + IT_0099*IT_0270 + IT_0221
      *IT_0271 + IT_0240*IT_0272 + IT_0250*IT_0273 + conj(IT_0096)*IT_0275);
    return create_ccomplex_return(IT_0276);
}

