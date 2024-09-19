#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_L_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_L_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
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
    const creal_t m_ss_L = param->m_ss_L;
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
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = pow(m_W, -2);
    const ccomplex_t IT_0010 = pow(m_c, 2);
    const ccomplex_t IT_0011 = pow(m_t, 2);
    const ccomplex_t IT_0012 = pow(m_u, 2);
    const ccomplex_t IT_0013 = pow(V_us, 2);
    const ccomplex_t IT_0014 = cpow(IT_0001, -2);
    const ccomplex_t IT_0015 = IT_0002*(IT_0004 + (-3)*IT_0006) + -IT_0008*
      (IT_0004 + (-3)*IT_0006*(1 + (-2)*IT_0009*(V_cs*conj(V_cs)*IT_0010 + V_ts
      *conj(V_ts)*IT_0011 + IT_0012*IT_0013)*IT_0014));
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0000*IT_0015;
    const ccomplex_t IT_0017 = 0.0833333333333333*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = pow(m_s, 2);
    const ccomplex_t IT_0020 = cpow(IT_0007, -1);
    const ccomplex_t IT_0021 = cpow(IT_0001, -1);
    const ccomplex_t IT_0022 = cpow(IT_0005, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*V_us*e_em*
      (IT_0001*IT_0018*IT_0019*IT_0020 + (-2)*IT_0007*(m_W*IT_0001 + (-0.5)
      *IT_0012*IT_0018*IT_0021))*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = cpow(IT_0024, 2);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = pow(m_Hp, 2);
    const ccomplex_t IT_0028 = pow(m_ss_L, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0027 + IT_0028 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *mu_h*IT_0018*IT_0022;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*mu_h*IT_0018*IT_0022;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0027 + IT_0028 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *mu_h*IT_0018*IT_0022;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = cpow(IT_0040, 2);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0027 + IT_0028 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = A_t*IT_0007*conj(U_st_10) + (-1.4142135623731)
      *e_em*IT_0022*(m_W*IT_0001*IT_0007*conj(U_st_00) + (-0.5)*IT_0018*(IT_0001
      *IT_0019*IT_0020*conj(U_st_00) + m_t*(m_t*IT_0007*IT_0021*conj(U_st_00) +
       mu_h*conj(U_st_10))));
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*V_ts*IT_0045;
    const ccomplex_t IT_0047 = A_t*IT_0007*U_st_10 + (-1.4142135623731)*e_em
      *IT_0022*(m_W*IT_0001*IT_0007*U_st_00 + (-0.5)*IT_0018*(IT_0001*IT_0019
      *IT_0020*U_st_00 + m_t*(m_t*IT_0007*IT_0021*U_st_00 + mu_h*U_st_10)));
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*conj(V_ts)*IT_0047;
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0027 + IT_0028 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = A_t*IT_0007*conj(U_st_11) + (-1.4142135623731)
      *e_em*IT_0022*(m_W*IT_0001*IT_0007*conj(U_st_01) + (-0.5)*IT_0018*(IT_0001
      *IT_0019*IT_0020*conj(U_st_01) + m_t*(m_t*IT_0007*IT_0021*conj(U_st_01) +
       mu_h*conj(U_st_11))));
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*V_ts*IT_0053;
    const ccomplex_t IT_0055 = A_t*IT_0007*U_st_11 + (-1.4142135623731)*e_em
      *IT_0022*(m_W*IT_0001*IT_0007*U_st_01 + (-0.5)*IT_0018*(IT_0001*IT_0019
      *IT_0020*U_st_01 + m_t*(m_t*IT_0007*IT_0021*U_st_01 + mu_h*U_st_11)));
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*conj(V_ts)*IT_0055;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0027 + IT_0028 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = sin(alpha);
    const ccomplex_t IT_0062 = cos(alpha);
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*e_em*(IT_0018*IT_0019
      *IT_0020*IT_0022*IT_0061 + (-0.166666666666667)*m_W*(IT_0004*IT_0005 + 3
      *IT_0022)*(IT_0007*IT_0061 + IT_0001*IT_0062));
    const ccomplex_t IT_0064 = IT_0004*IT_0005;
    const ccomplex_t IT_0065 = -IT_0022 + IT_0064;
    const ccomplex_t IT_0066 = cpow(IT_0001, 3);
    const ccomplex_t IT_0067 = 3*IT_0022 + IT_0064;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*m_W*e_em*(IT_0062*IT_0065
      *IT_0066 + IT_0007*(IT_0002*IT_0061*IT_0067 + -IT_0007*(IT_0007*IT_0061
      *IT_0065 + IT_0001*IT_0062*IT_0067)));
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = IT_0063*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0028 + -reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*e_em*(IT_0018*IT_0019
      *IT_0020*IT_0022*IT_0062 + 0.166666666666667*m_W*(IT_0004*IT_0005 + 3
      *IT_0022)*(IT_0001*IT_0061 + -IT_0007*IT_0062));
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = -IT_0022;
    const ccomplex_t IT_0077 = IT_0064 + IT_0076;
    const ccomplex_t IT_0078 = -IT_0007;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*m_W*e_em*(IT_0061*IT_0066
      *IT_0077 + IT_0078*(IT_0002*IT_0062*IT_0067 + -(IT_0001*IT_0061*IT_0067 + 
      -IT_0007*IT_0062*IT_0077)*IT_0078));
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = IT_0075*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_12 + (-2)*IT_0028 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*(IT_0001*IT_0018*IT_0019*IT_0020 + (-2)*IT_0007*(m_W*IT_0001 + (-0.5
      )*IT_0010*IT_0018*IT_0021))*IT_0022;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em*
      (IT_0001*IT_0018*IT_0019*IT_0020 + (-2)*IT_0007*(m_W*IT_0001 + (-0.5)
      *IT_0010*IT_0018*IT_0021))*IT_0022;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = IT_0086*IT_0088;
    const ccomplex_t IT_0090 = cpow((-2)*s_13 + IT_0027 + IT_0028 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0091 = IT_0089*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0094 = 0.666666666666667*IT_0093;
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = s_12 + IT_0028;
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = 2*IT_0093;
    const ccomplex_t IT_0099 = 0.333333333333333*IT_0093;
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = s_13 + s_23;
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = IT_0097 + IT_0102;
    const ccomplex_t IT_0104 = IT_0094*IT_0098;
    const ccomplex_t IT_0105 = s_23*IT_0104;
    const ccomplex_t IT_0106 = IT_0093*IT_0099;
    const ccomplex_t IT_0107 = IT_0096*IT_0106;
    const ccomplex_t IT_0108 = -IT_0105 + (-2)*IT_0107;
    const ccomplex_t IT_0109 = IT_0103 + IT_0108;
    const ccomplex_t IT_0110 = cpow(s_12 + IT_0028 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = cpow(IT_0003, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*-1)*e_em*(IT_0003*IT_0022 + 
      -IT_0005*IT_0113);
    const ccomplex_t IT_0115 = (0 + _Complex_I*3)*e_em*(IT_0003*IT_0022 + 1./3
      *IT_0005*IT_0113);
    const ccomplex_t IT_0116 = 0.166666666666667*IT_0115;
    const ccomplex_t IT_0117 = IT_0114*IT_0116;
    const ccomplex_t IT_0118 = s_14 + s_24;
    const ccomplex_t IT_0119 = IT_0117*IT_0118;
    const ccomplex_t IT_0120 = 0.5*IT_0114;
    const ccomplex_t IT_0121 = 0.333333333333333*IT_0115;
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = IT_0096*IT_0122;
    const ccomplex_t IT_0124 = IT_0119 + IT_0123;
    const ccomplex_t IT_0125 = IT_0116*IT_0120;
    const ccomplex_t IT_0126 = IT_0096*IT_0125;
    const ccomplex_t IT_0127 = IT_0114*IT_0121;
    const ccomplex_t IT_0128 = s_24*IT_0127;
    const ccomplex_t IT_0129 = (-2)*IT_0126 + -IT_0128;
    const ccomplex_t IT_0130 = IT_0124 + IT_0129;
    const ccomplex_t IT_0131 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0028 + -reg_prop, -1);
    const ccomplex_t IT_0132 = IT_0130*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = IT_0017 + IT_0030 + IT_0038 + IT_0044 + IT_0052
       + IT_0060 + -IT_0073 + -IT_0084 + IT_0092 + (-0.5)*IT_0112 + IT_0133;
    return create_ccomplex_return(3*IT_0134*conj(IT_0134));
}

