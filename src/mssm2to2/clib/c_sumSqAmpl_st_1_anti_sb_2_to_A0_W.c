#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_sb_2_to_A0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_sb_2_to_A0_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0001*U_st_00;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = 1.4142135623731*A_b*V_tb*V_ts*IT_0004*conj
      (U_sb_11);
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = pow(m_A0, 2);
    const ccomplex_t IT_0008 = pow(m_sb_2, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0003*IT_0011;
    const ccomplex_t IT_0013 = e_em*IT_0001;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = cpow(IT_0004, -1);
    const ccomplex_t IT_0018 = cpow(IT_0015, -1);
    const ccomplex_t IT_0019 = m_t*conj(U_sb_01)*(m_t*IT_0015*IT_0017*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0004*IT_0018*conj(U_sb_01)*U_st_00 + conj
      (U_sb_11)*(mu_h*U_st_00 + m_t*IT_0017*IT_0018*U_st_10));
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0004*conj
      (U_sb_11)*U_st_00 + (-1.4142135623731)*e_em*IT_0001*((-0.5)*IT_0016
      *IT_0019 + m_W*IT_0004*IT_0015*conj(U_sb_01)*U_st_00) + A_t*IT_0015*conj
      (U_sb_01)*U_st_10);
    const ccomplex_t IT_0021 = pow(m_st_1, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + -IT_0008 + -IT_0021 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0014*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_01)*U_st_00;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = m_b*e_em*mu_h*IT_0016*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0029 = IT_0001*IT_0028;
    const ccomplex_t IT_0030 = m_b*e_em*mu_h*IT_0016*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0031 = IT_0001*IT_0030;
    const ccomplex_t IT_0032 = cpow(V_tb, 2);
    const ccomplex_t IT_0033 = A_b*IT_0032*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0034 = IT_0004*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = A_b*IT_0032*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0037 = IT_0004*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (-0.5)*IT_0029 + 0.5*IT_0031 + (-0.5)*IT_0035 +
       0.5*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0007 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0027*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = m_b*e_em*mu_h*IT_0016*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0047 = IT_0001*IT_0046;
    const ccomplex_t IT_0048 = m_b*e_em*mu_h*IT_0016*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0049 = IT_0001*IT_0048;
    const ccomplex_t IT_0050 = A_b*IT_0032*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = A_b*IT_0032*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0054 = IT_0004*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (-0.5)*IT_0047 + 0.5*IT_0049 + (-0.5)*IT_0052 +
       0.5*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0045*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0001*U_st_00;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*A_b*V_tb*V_td*IT_0004*conj
      (U_sb_11);
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0062*IT_0067;
    const ccomplex_t IT_0069 = -IT_0012 + IT_0025 + -IT_0043 + -IT_0060 + 
      -IT_0068;
    const ccomplex_t IT_0070 = pow(m_W, -2);
    const ccomplex_t IT_0071 = pow(s_34, 2);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = m_t*e_em*mu_h*IT_0016*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0074 = IT_0001*IT_0073;
    const ccomplex_t IT_0075 = m_t*e_em*mu_h*IT_0016*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0076 = IT_0001*IT_0075;
    const ccomplex_t IT_0077 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0078 = IT_0015*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0081 = IT_0015*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (-0.5)*IT_0074 + 0.5*IT_0076 + (-0.5)*IT_0079 +
       0.5*IT_0082;
    const ccomplex_t IT_0084 = cpow((-2)*s_13 + IT_0007 + (0 + _Complex_I*1)
      *m_st_1*Gamma_t1 + reg_prop, -1);
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0027*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_01)*U_st_01;
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = m_t*e_em*mu_h*IT_0016*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0091 = IT_0001*IT_0090;
    const ccomplex_t IT_0092 = m_t*e_em*mu_h*IT_0016*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0093 = IT_0001*IT_0092;
    const ccomplex_t IT_0094 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0095 = IT_0015*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0098 = IT_0015*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (-0.5)*IT_0091 + 0.5*IT_0093 + (-0.5)*IT_0096 +
       0.5*IT_0099;
    const ccomplex_t IT_0101 = cpow((-2)*s_13 + IT_0007 + IT_0021 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0089*IT_0103;
    const ccomplex_t IT_0105 = IT_0012 + IT_0043 + IT_0060 + IT_0068 + IT_0087
       + IT_0104;
    const ccomplex_t IT_0106 = s_24*s_34*IT_0070;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = s_23 + IT_0107;
    const ccomplex_t IT_0109 = pow(s_24, 2);
    const ccomplex_t IT_0110 = IT_0070*IT_0109;
    return create_ccomplex_return((-3)*IT_0069*(conj(IT_0069)*(IT_0007 + 
      -IT_0072) + conj(IT_0105)*IT_0108) + (-3)*IT_0105*(conj(IT_0069)*IT_0108 +
       conj(IT_0105)*(IT_0008 + -IT_0110)));
}

