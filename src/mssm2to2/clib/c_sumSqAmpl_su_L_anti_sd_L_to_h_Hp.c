#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sd_L_to_h_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sd_L_to_h_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = pow(m_W, -2);
    const ccomplex_t IT_0004 = pow(m_u, 2);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = pow(m_d, 2);
    const ccomplex_t IT_0008 = cpow(IT_0001, -2);
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*V_ud*IT_0000*
      (IT_0001*(IT_0002 + -IT_0002*IT_0003*IT_0004*IT_0006) + -IT_0005*(1 + 
      -IT_0003*IT_0007*IT_0008)*IT_0009)*IT_0011;
    const ccomplex_t IT_0013 = (-0.25)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = cpow(IT_0001, -1);
    const ccomplex_t IT_0016 = cpow(IT_0005, -1);
    const ccomplex_t IT_0017 = cpow(IT_0010, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em*
      (IT_0005*IT_0007*IT_0014*IT_0015 + (-2)*IT_0001*(m_W*IT_0005 + (-0.5)
      *IT_0004*IT_0014*IT_0016))*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = cos(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -2);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0004
      *IT_0014*IT_0016*IT_0017 + (-0.5)*m_W*(IT_0002*IT_0005 + IT_0001*IT_0009)*
      (IT_0017 + -1./3*IT_0010*IT_0021));
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = IT_0019*IT_0023;
    const ccomplex_t IT_0025 = pow(m_h, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_13 + IT_0025 + (0 + _Complex_I*1)
      *m_su_L*Gamma_ul + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*e_em*(IT_0007*IT_0009
      *IT_0014*IT_0015*IT_0017 + (-0.5)*m_W*(IT_0002*IT_0005 + IT_0001*IT_0009)*
      (IT_0017 + 1./3*IT_0010*IT_0021));
    const ccomplex_t IT_0030 = IT_0019*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0025 + (0 + _Complex_I*1)
      *m_sd_L*Gamma_dl + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = cpow(IT_0005, 3);
    const ccomplex_t IT_0035 = cpow(IT_0001, 2);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*m_W*e_em*(IT_0009*IT_0017
      *IT_0034 + (-2)*IT_0001*(IT_0005*(IT_0002*IT_0005 + IT_0001*IT_0009)*(0.5
      *IT_0017 + IT_0010*IT_0021) + -1./2*IT_0002*IT_0017*IT_0035));
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = cpow(IT_0005, 2);
    const ccomplex_t IT_0039 = (0 + _Complex_I*(-1.4142135623731))*V_ud*e_em
      *IT_0017*((IT_0004 + -IT_0007)*IT_0014 + m_W*(IT_0035 + -IT_0038));
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = IT_0037*IT_0040;
    const ccomplex_t IT_0042 = pow(m_W, 2);
    const ccomplex_t IT_0043 = pow(m_sd_L, 2);
    const ccomplex_t IT_0044 = pow(m_su_L, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_12 + IT_0042 + -IT_0043 + -IT_0044 
      + -reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0041*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0009
      *IT_0014*IT_0016*IT_0017;
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *mu_h*IT_0014*IT_0017;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_13 + IT_0025 + IT_0044 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0014*IT_0017;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0009 + m_d*e_em*mu_h*IT_0002*IT_0014*IT_0015*IT_0017);
    const ccomplex_t IT_0059 = 0.5*IT_0058;
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_23 + IT_0025 + IT_0043 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *mu_h*IT_0014*IT_0017;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*A_b*conj(V_td)
      *V_ts*IT_0009;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0025 + IT_0043 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0017
      *IT_0072*(m_W*IT_0001*IT_0005*conj(U_sb_01) + (-0.5)*IT_0014*(IT_0001
      *IT_0004*IT_0016*conj(U_sb_01) + m_b*(m_b*IT_0005*IT_0015*conj(U_sb_01) +
       mu_h*conj(U_sb_11))))*V_ub_mod;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0009*U_sb_11;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0025 + IT_0043 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0025 + IT_0043 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0017
      *IT_0072*(m_W*IT_0001*IT_0005*conj(U_sb_00) + (-0.5)*IT_0014*(IT_0001
      *IT_0004*IT_0016*conj(U_sb_00) + m_b*(m_b*IT_0005*IT_0015*conj(U_sb_00) +
       mu_h*conj(U_sb_10))))*V_ub_mod;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0009*U_sb_10;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = IT_0081*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0043 + -IT_0044 + -reg_prop, -1);
    const ccomplex_t IT_0090 = s_12 + IT_0043;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0017;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 +
       IT_0005*IT_0009)*IT_0017;
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = IT_0092*IT_0094;
    const ccomplex_t IT_0096 = IT_0090*IT_0095;
    const ccomplex_t IT_0097 = -IT_0093;
    const ccomplex_t IT_0098 = IT_0092*IT_0097;
    const ccomplex_t IT_0099 = -s_23;
    const ccomplex_t IT_0100 = IT_0090 + IT_0099;
    const ccomplex_t IT_0101 = IT_0098*IT_0100;
    const ccomplex_t IT_0102 = (-0.5)*IT_0091;
    const ccomplex_t IT_0103 = IT_0097*IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*s_13 + (-0.5)*s_23 + 0.5*IT_0043 + 0.5
      *IT_0044;
    const ccomplex_t IT_0105 = s_12 + IT_0104;
    const ccomplex_t IT_0106 = IT_0103*IT_0105;
    const ccomplex_t IT_0107 = (-0.5)*IT_0096 + 0.5*IT_0101 + -IT_0106;
    const ccomplex_t IT_0108 = IT_0094*IT_0102;
    const ccomplex_t IT_0109 = 0.5*IT_0043 + 0.5*IT_0044;
    const ccomplex_t IT_0110 = s_12 + IT_0109;
    const ccomplex_t IT_0111 = IT_0108*IT_0110;
    const ccomplex_t IT_0112 = IT_0107 + IT_0111;
    const ccomplex_t IT_0113 = IT_0089*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0010*IT_0021;
    const ccomplex_t IT_0116 = -IT_0017 + IT_0115;
    const ccomplex_t IT_0117 = 3*IT_0017 + IT_0115;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_W*e_em*(IT_0002*IT_0034
      *IT_0116 + IT_0001*(IT_0009*IT_0038*IT_0117 + -IT_0001*(IT_0001*IT_0009
      *IT_0116 + IT_0002*IT_0005*IT_0117)));
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0019*IT_0119;
    const ccomplex_t IT_0121 = cpow((-2)*s_12 + -IT_0043 + -IT_0044 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = IT_0013 + IT_0028 + IT_0033 + -IT_0047 +
       IT_0055 + IT_0063 + IT_0071 + IT_0080 + IT_0088 + (-2)*IT_0114 + -IT_0123;
    return create_ccomplex_return(3*IT_0124*conj(IT_0124));
}

