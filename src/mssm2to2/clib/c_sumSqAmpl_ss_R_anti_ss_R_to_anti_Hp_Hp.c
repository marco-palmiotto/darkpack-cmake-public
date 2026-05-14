#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_anti_ss_R_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_anti_ss_R_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_R = param->m_ss_R;
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
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_ts = param->V_ts;
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
    const ccomplex_t IT_0008 = pow(m_s, 2);
    const ccomplex_t IT_0009 = cpow(IT_0003, -2);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0000*(IT_0002*(IT_0004 + 
      -IT_0006) + 3*IT_0006*IT_0007*IT_0008*IT_0009*IT_0011);
    const ccomplex_t IT_0013 = (-0.166666666666667)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = cpow(IT_0003, -1);
    const ccomplex_t IT_0016 = cpow(IT_0005, -1);
    const ccomplex_t IT_0017 = cpow(IT_0010, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*m_s*m_u*V_us
      *e_em*IT_0014*IT_0015*IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0019, 2);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = pow(m_Hp, 2);
    const ccomplex_t IT_0023 = pow(m_ss_R, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*V_ts*(A_b*IT_0005*conj
      (U_st_00) + 0.707106781186548*m_s*e_em*IT_0014*IT_0017*(mu_h*conj(U_st_00)
       + m_t*IT_0015*IT_0016*conj(U_st_10)));
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*conj(V_ts)*(A_b*IT_0005
      *U_st_00 + 0.707106781186548*m_s*e_em*IT_0014*IT_0017*(mu_h*U_st_00 + m_t
      *IT_0015*IT_0016*U_st_10));
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*V_ts*(A_b*IT_0005*conj
      (U_st_01) + 0.707106781186548*m_s*e_em*IT_0014*IT_0017*(mu_h*conj(U_st_01)
       + m_t*IT_0015*IT_0016*conj(U_st_11)));
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*conj(V_ts)*(A_b*IT_0005
      *U_st_01 + 0.707106781186548*m_s*e_em*IT_0014*IT_0017*(mu_h*U_st_01 + m_t
      *IT_0015*IT_0016*U_st_11));
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow(s_12 + IT_0023 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0040 = 0.666666666666667*IT_0039;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = s_12 + IT_0023;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = s_13 + s_23;
    const ccomplex_t IT_0045 = 0.333333333333333*IT_0039;
    const ccomplex_t IT_0046 = 2*IT_0039;
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = IT_0044*IT_0047;
    const ccomplex_t IT_0049 = IT_0043 + IT_0048;
    const ccomplex_t IT_0050 = IT_0039*IT_0045;
    const ccomplex_t IT_0051 = IT_0042*IT_0050;
    const ccomplex_t IT_0052 = IT_0040*IT_0046;
    const ccomplex_t IT_0053 = s_23*IT_0052;
    const ccomplex_t IT_0054 = (-2)*IT_0051 + -IT_0053;
    const ccomplex_t IT_0055 = IT_0049 + IT_0054;
    const ccomplex_t IT_0056 = IT_0038*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = cpow(IT_0001, -1);
    const ccomplex_t IT_0059 = (0 + _Complex_I*-1)*e_em*(IT_0001*IT_0017 + 
      -IT_0010*IT_0058);
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*e_em*IT_0010*IT_0058;
    const ccomplex_t IT_0061 = (-0.666666666666667)*IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = s_24*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0059;
    const ccomplex_t IT_0065 = IT_0061*IT_0064;
    const ccomplex_t IT_0066 = IT_0042*IT_0065;
    const ccomplex_t IT_0067 = (-0.333333333333333)*IT_0060;
    const ccomplex_t IT_0068 = IT_0059*IT_0067;
    const ccomplex_t IT_0069 = s_14 + s_24;
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = IT_0064*IT_0067;
    const ccomplex_t IT_0072 = IT_0042*IT_0071;
    const ccomplex_t IT_0073 = -IT_0066 + -IT_0070 + 2*IT_0072;
    const ccomplex_t IT_0074 = IT_0063 + IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = cos(alpha);
    const ccomplex_t IT_0079 = sin(alpha);
    const ccomplex_t IT_0080 = (0 + _Complex_I*-3)*e_em*(IT_0008*IT_0014
      *IT_0015*IT_0017*IT_0078 + -1./3*m_W*IT_0002*IT_0010*(IT_0003*IT_0078 + 
      -IT_0005*IT_0079));
    const ccomplex_t IT_0081 = 0.333333333333333*IT_0080;
    const ccomplex_t IT_0082 = -IT_0017;
    const ccomplex_t IT_0083 = IT_0002*IT_0010;
    const ccomplex_t IT_0084 = IT_0082 + IT_0083;
    const ccomplex_t IT_0085 = cpow(IT_0005, 3);
    const ccomplex_t IT_0086 = -IT_0003;
    const ccomplex_t IT_0087 = 3*IT_0017 + IT_0083;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*m_W*e_em*(IT_0079*IT_0084
      *IT_0085 + IT_0086*(IT_0006*IT_0078*IT_0087 + IT_0086*(IT_0003*IT_0078
      *IT_0084 + -IT_0005*IT_0079*IT_0087)));
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = IT_0081*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_12 + (-2)*IT_0023 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0014*IT_0017;
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *mu_h*IT_0014*IT_0017;
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = IT_0095*IT_0097;
    const ccomplex_t IT_0099 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *mu_h*IT_0014*IT_0017;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = cpow(IT_0103, 2);
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0107 = IT_0105*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*-3)*e_em*(IT_0008*IT_0014
      *IT_0015*IT_0017*IT_0079 + -1./3*m_W*IT_0002*IT_0010*(IT_0005*IT_0078 +
       IT_0003*IT_0079));
    const ccomplex_t IT_0109 = (-0.333333333333333)*IT_0108;
    const ccomplex_t IT_0110 = -IT_0017 + IT_0083;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*m_W*e_em*(IT_0078*IT_0085
      *IT_0110 + IT_0003*(IT_0006*IT_0079*IT_0087 + -IT_0003*(IT_0005*IT_0078
      *IT_0087 + IT_0003*IT_0079*IT_0110)));
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = IT_0109*IT_0112;
    const ccomplex_t IT_0114 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*m_c*m_s*V_cs
      *e_em*IT_0014*IT_0015*IT_0016*IT_0017;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_c*m_s*conj
      (V_cs)*e_em*IT_0014*IT_0015*IT_0016*IT_0017;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = IT_0118*IT_0120;
    const ccomplex_t IT_0122 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0123 = IT_0121*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = IT_0013 + IT_0025 + IT_0031 + IT_0037 + (-0.5)
      *IT_0057 + -IT_0077 + -IT_0093 + IT_0101 + IT_0107 + -IT_0116 + IT_0124;
    return create_ccomplex_return(3*IT_0125*conj(IT_0125));
}

