#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_anti_sb_2_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_anti_sb_2_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = pow(m_W, -2);
    const ccomplex_t IT_0008 = pow(m_b, 2);
    const ccomplex_t IT_0009 = cpow(IT_0003, -2);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = pow(m_c, 2);
    const ccomplex_t IT_0013 = pow(V_cb, 2);
    const ccomplex_t IT_0014 = pow(m_t, 2);
    const ccomplex_t IT_0015 = cpow(V_tb, 2);
    const ccomplex_t IT_0016 = pow(m_u, 2);
    const ccomplex_t IT_0017 = pow(V_ub_mod, 2);
    const ccomplex_t IT_0018 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0019 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0020 = cpow(IT_0005, -2);
    const ccomplex_t IT_0021 = IT_0006*(IT_0002 + (-3)*IT_0011) + -IT_0004*
      (IT_0002 + (-3)*IT_0011*(1 + (-2)*IT_0007*(IT_0012*IT_0013 + IT_0014
      *IT_0015 + IT_0016*IT_0017*IT_0018*IT_0019)*IT_0020));
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0000*((-0.5)*IT_0021
      *U_sb_01*conj(U_sb_01) + (IT_0002*(IT_0004 + -IT_0006) + 3*IT_0006*IT_0007
      *IT_0008*IT_0009*IT_0011)*U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0023 = (-0.166666666666667)*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0025 = 0.666666666666667*IT_0024;
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = pow(m_sb_2, 2);
    const ccomplex_t IT_0028 = s_12 + IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = 2*IT_0024;
    const ccomplex_t IT_0031 = 0.333333333333333*IT_0024;
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = s_13 + s_23;
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = IT_0029 + IT_0034;
    const ccomplex_t IT_0036 = IT_0025*IT_0030;
    const ccomplex_t IT_0037 = s_23*IT_0036;
    const ccomplex_t IT_0038 = IT_0024*IT_0031;
    const ccomplex_t IT_0039 = IT_0028*IT_0038;
    const ccomplex_t IT_0040 = -IT_0037 + (-2)*IT_0039;
    const ccomplex_t IT_0041 = IT_0035 + IT_0040;
    const ccomplex_t IT_0042 = cpow(s_12 + IT_0027 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0027 + -reg_prop, -1);
    const ccomplex_t IT_0046 = cpow(IT_0001, -1);
    const ccomplex_t IT_0047 = cpow(IT_0010, -1);
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*e_em*(IT_0010*IT_0046 + 
      -IT_0001*IT_0047);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0010
      *IT_0046 + 3*IT_0001*IT_0047)*U_sb_01*conj(U_sb_01) + IT_0010*IT_0046
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0050 = (-0.333333333333333)*IT_0049;
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*s_13 + (-0.5)*s_23;
    const ccomplex_t IT_0053 = IT_0028 + IT_0052;
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = 0.5*IT_0048;
    const ccomplex_t IT_0056 = IT_0050*IT_0055;
    const ccomplex_t IT_0057 = IT_0028*IT_0056;
    const ccomplex_t IT_0058 = (-0.666666666666667)*IT_0049;
    const ccomplex_t IT_0059 = IT_0055*IT_0058;
    const ccomplex_t IT_0060 = IT_0028*IT_0059;
    const ccomplex_t IT_0061 = IT_0048*IT_0058;
    const ccomplex_t IT_0062 = -s_23;
    const ccomplex_t IT_0063 = IT_0028 + IT_0062;
    const ccomplex_t IT_0064 = IT_0061*IT_0063;
    const ccomplex_t IT_0065 = -IT_0057 + 0.5*IT_0060 + (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0054 + IT_0065;
    const ccomplex_t IT_0067 = IT_0045*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0002*IT_0010;
    const ccomplex_t IT_0070 = -IT_0047 + IT_0069;
    const ccomplex_t IT_0071 = cos(alpha);
    const ccomplex_t IT_0072 = cpow(IT_0005, 3);
    const ccomplex_t IT_0073 = 3*IT_0047 + IT_0069;
    const ccomplex_t IT_0074 = sin(alpha);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_W*e_em*(IT_0070*IT_0071
      *IT_0072 + IT_0003*(IT_0006*IT_0073*IT_0074 + -IT_0003*(IT_0005*IT_0071
      *IT_0073 + IT_0003*IT_0070*IT_0074)));
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = pow(m_W, -1);
    const ccomplex_t IT_0078 = m_b*mu_h*IT_0077;
    const ccomplex_t IT_0079 = cpow(IT_0003, -1);
    const ccomplex_t IT_0080 = IT_0003*IT_0074;
    const ccomplex_t IT_0081 = IT_0005*IT_0071;
    const ccomplex_t IT_0082 = IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = 2*U_sb_01*((IT_0008*IT_0047*IT_0074*IT_0077
      *IT_0079 + (-0.166666666666667)*m_W*(IT_0002*IT_0010 + 3*IT_0047)*IT_0082)
      *conj(U_sb_01) + 1./2*IT_0047*IT_0071*IT_0078*IT_0079*conj(U_sb_11)) +
       U_sb_11*(IT_0047*IT_0071*IT_0078*IT_0079*conj(U_sb_01) + 2*(IT_0008
      *IT_0047*IT_0074*IT_0077*IT_0079 + -1./3*m_W*IT_0002*IT_0010*IT_0082)*conj
      (U_sb_11));
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(e_em*IT_0083 +
       1.4142135623731*A_b*IT_0015*IT_0074*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0076*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0027 + -reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = -IT_0047;
    const ccomplex_t IT_0091 = IT_0069 + IT_0090;
    const ccomplex_t IT_0092 = -IT_0003;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_W*e_em*(IT_0072*IT_0074
      *IT_0091 + IT_0092*(IT_0006*IT_0071*IT_0073 + -(IT_0005*IT_0073*IT_0074 + 
      -IT_0003*IT_0071*IT_0091)*IT_0092));
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = -IT_0005*IT_0074;
    const ccomplex_t IT_0096 = IT_0003*IT_0071;
    const ccomplex_t IT_0097 = IT_0095 + IT_0096;
    const ccomplex_t IT_0098 = (-2)*U_sb_01*((IT_0008*IT_0047*IT_0071*IT_0077
      *IT_0079 + (-0.166666666666667)*m_W*(IT_0002*IT_0010 + 3*IT_0047)*IT_0097)
      *conj(U_sb_01) + -1./2*IT_0047*IT_0074*IT_0078*IT_0079*conj(U_sb_11)) +
       U_sb_11*(IT_0047*IT_0074*IT_0078*IT_0079*conj(U_sb_01) + (-2)*(IT_0008
      *IT_0047*IT_0071*IT_0077*IT_0079 + -1./3*m_W*IT_0002*IT_0010*IT_0097)*conj
      (U_sb_11));
    const ccomplex_t IT_0099 = (0 + _Complex_I*-1)*(e_em*IT_0098 + (
      -1.4142135623731)*A_b*IT_0015*IT_0071*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0094*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_12 + (-2)*IT_0027 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = cpow(IT_0005, -1);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0047*(m_W*IT_0003*IT_0005*conj(U_sb_01) + (-0.5)*IT_0077*(IT_0003
      *IT_0012*IT_0105*conj(U_sb_01) + m_b*(m_b*IT_0005*IT_0079*conj(U_sb_01) +
       mu_h*conj(U_sb_11))));
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0047*(m_W*IT_0003*IT_0005*U_sb_01 + (-0.5)*IT_0077*(IT_0003*IT_0012
      *IT_0105*U_sb_01 + m_b*(m_b*IT_0005*IT_0079*U_sb_01 + mu_h*U_sb_11)));
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = IT_0107*IT_0109;
    const ccomplex_t IT_0111 = pow(m_Hp, 2);
    const ccomplex_t IT_0112 = cpow((-2)*s_13 + IT_0027 + IT_0111 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0113 = IT_0110*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0018
      *IT_0047*(m_W*IT_0003*IT_0005*conj(U_sb_01) + (-0.5)*IT_0077*(IT_0003
      *IT_0016*IT_0105*conj(U_sb_01) + m_b*(m_b*IT_0005*IT_0079*conj(U_sb_01) +
       mu_h*conj(U_sb_11))))*V_ub_mod;
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0019
      *IT_0047*(m_W*IT_0003*IT_0005*U_sb_01 + (-0.5)*IT_0077*(IT_0003*IT_0016
      *IT_0105*U_sb_01 + m_b*(m_b*IT_0005*IT_0079*U_sb_01 + mu_h*U_sb_11)))
      *V_ub_mod;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = IT_0116*IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_13 + IT_0027 + IT_0111 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0121 = IT_0119*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0018*IT_0047*IT_0077*(mu_h*conj(U_sb_01) + m_b*IT_0079*IT_0105*conj
      (U_sb_11))*V_ub_mod;
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0019*IT_0047*IT_0077*(mu_h*U_sb_01 + m_b*IT_0079*IT_0105*U_sb_11)
      *V_ub_mod;
    const ccomplex_t IT_0126 = 0.5*IT_0125;
    const ccomplex_t IT_0127 = IT_0124*IT_0126;
    const ccomplex_t IT_0128 = cpow((-2)*s_13 + IT_0027 + IT_0111 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0129 = IT_0127*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = m_t*conj(U_sb_01)*(m_t*IT_0003*IT_0105*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0005*IT_0079*conj(U_sb_01)*U_st_00 + conj
      (U_sb_11)*(mu_h*U_st_00 + m_t*IT_0079*IT_0105*U_st_10));
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_11)*U_st_00 + (-1.4142135623731)*e_em*IT_0047*((-0.5)*IT_0077
      *IT_0131 + m_W*IT_0003*IT_0005*conj(U_sb_01)*U_st_00) + A_t*IT_0003*conj
      (U_sb_01)*U_st_10);
    const ccomplex_t IT_0133 = m_t*U_sb_01*(m_t*IT_0003*IT_0105*conj(U_st_00) 
      + mu_h*conj(U_st_10)) + m_b*(m_b*IT_0005*IT_0079*U_sb_01*conj(U_st_00) +
       U_sb_11*(mu_h*conj(U_st_00) + m_t*IT_0079*IT_0105*conj(U_st_10)));
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*U_sb_11
      *conj(U_st_00) + (-1.4142135623731)*e_em*IT_0047*((-0.5)*IT_0077*IT_0133 +
       m_W*IT_0003*IT_0005*U_sb_01*conj(U_st_00)) + A_t*IT_0003*U_sb_01*conj
      (U_st_10));
    const ccomplex_t IT_0135 = IT_0132*IT_0134;
    const ccomplex_t IT_0136 = cpow((-2)*s_13 + IT_0027 + IT_0111 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0137 = IT_0135*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = m_t*conj(U_sb_01)*(m_t*IT_0003*IT_0105*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0005*IT_0079*conj(U_sb_01)*U_st_01 + conj
      (U_sb_11)*(mu_h*U_st_01 + m_t*IT_0079*IT_0105*U_st_11));
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_11)*U_st_01 + (-1.4142135623731)*e_em*IT_0047*((-0.5)*IT_0077
      *IT_0139 + m_W*IT_0003*IT_0005*conj(U_sb_01)*U_st_01) + A_t*IT_0003*conj
      (U_sb_01)*U_st_11);
    const ccomplex_t IT_0141 = m_t*U_sb_01*(m_t*IT_0003*IT_0105*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0005*IT_0079*U_sb_01*conj(U_st_01) +
       U_sb_11*(mu_h*conj(U_st_01) + m_t*IT_0079*IT_0105*conj(U_st_11)));
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*U_sb_11
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0047*((-0.5)*IT_0077*IT_0141 +
       m_W*IT_0003*IT_0005*U_sb_01*conj(U_st_01)) + A_t*IT_0003*U_sb_01*conj
      (U_st_11));
    const ccomplex_t IT_0143 = IT_0140*IT_0142;
    const ccomplex_t IT_0144 = cpow((-2)*s_13 + IT_0027 + IT_0111 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0145 = IT_0143*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0047*IT_0077*(mu_h*conj(U_sb_01) + m_b*IT_0079*IT_0105*conj(U_sb_11));
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0047*IT_0077*(mu_h*U_sb_01 + m_b*IT_0079*IT_0105*U_sb_11);
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = IT_0148*IT_0150;
    const ccomplex_t IT_0152 = cpow((-2)*s_13 + IT_0027 + IT_0111 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0023 + (-0.5)*IT_0044 + 2*IT_0068 + -IT_0089
       + -IT_0104 + IT_0114 + IT_0122 + IT_0130 + IT_0138 + IT_0146 + IT_0154;
    return create_ccomplex_return(3*IT_0155*conj(IT_0155));
}

