#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_anti_Hp_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
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
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 +
       IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = m_t*mu_h*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0002, -1);
    const ccomplex_t IT_0011 = IT_0000*IT_0003;
    const ccomplex_t IT_0012 = IT_0001*IT_0002;
    const ccomplex_t IT_0013 = IT_0011 + IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -2);
    const ccomplex_t IT_0016 = pow(m_t, 2);
    const ccomplex_t IT_0017 = 2*U_st_01*(((-0.5)*m_W*IT_0013*(IT_0005 + -1./3
      *IT_0004*IT_0015) + IT_0001*IT_0005*IT_0008*IT_0010*IT_0016)*conj(U_st_01)
       + 1./2*IT_0003*IT_0005*IT_0009*IT_0010*conj(U_st_11)) + U_st_11*(IT_0003
      *IT_0005*IT_0009*IT_0010*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0004
      *IT_0013*IT_0015 + (-1.5)*IT_0001*IT_0005*IT_0008*IT_0010*IT_0016)*conj
      (U_st_11));
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(e_em*IT_0017 +
       1.4142135623731*A_t*IT_0001*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_st_2, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0020 + -reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0007*IT_0023;
    const ccomplex_t IT_0025 = e_em*IT_0005;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = m_t*e_em*mu_h*IT_0008*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0028 = IT_0005*IT_0027;
    const ccomplex_t IT_0029 = m_t*e_em*mu_h*IT_0008*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0030 = IT_0005*IT_0029;
    const ccomplex_t IT_0031 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0032 = IT_0000*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0035 = IT_0000*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0028 + 0.5*IT_0030 + (-0.5)*IT_0033 +
       0.5*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_12 + (-2)*IT_0020 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0026*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003)*IT_0005;
    const ccomplex_t IT_0043 = IT_0000*IT_0001;
    const ccomplex_t IT_0044 = -IT_0002*IT_0003 + IT_0043;
    const ccomplex_t IT_0045 = (-2)*U_st_01*((IT_0003*IT_0005*IT_0008*IT_0010
      *IT_0016 + 0.5*m_W*(IT_0005 + -1./3*IT_0004*IT_0015)*IT_0044)*conj(U_st_01
      ) + -1./2*IT_0001*IT_0005*IT_0009*IT_0010*conj(U_st_11)) + U_st_11*
      (IT_0001*IT_0005*IT_0009*IT_0010*conj(U_st_01) + (-1.33333333333333)*(1.5
      *IT_0003*IT_0005*IT_0008*IT_0010*IT_0016 + m_W*IT_0004*IT_0015*IT_0044)
      *conj(U_st_11));
    const ccomplex_t IT_0046 = (0 + _Complex_I*-1)*(e_em*IT_0045 + (
      -1.4142135623731)*A_t*IT_0003*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + (-2)*IT_0020 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0042*IT_0050;
    const ccomplex_t IT_0052 = -IT_0024 + -IT_0041 + -IT_0051;
    const ccomplex_t IT_0053 = pow(m_W, -2);
    const ccomplex_t IT_0054 = pow(s_14, 2);
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0005*U_st_01;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = pow(m_s, 2);
    const ccomplex_t IT_0059 = cpow(IT_0000, -1);
    const ccomplex_t IT_0060 = A_t*IT_0000*conj(U_st_11) + (-1.4142135623731)
      *e_em*IT_0005*(m_W*IT_0000*IT_0002*conj(U_st_01) + (-0.5)*IT_0008*(IT_0002
      *IT_0058*IT_0059*conj(U_st_01) + m_t*(m_t*IT_0000*IT_0010*conj(U_st_01) +
       mu_h*conj(U_st_11))));
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*V_ts*IT_0060;
    const ccomplex_t IT_0062 = pow(m_Hp, 2);
    const ccomplex_t IT_0063 = cpow((-2)*s_23 + IT_0020 + IT_0062 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0061*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0057*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0005*conj(U_sb_01)*U_st_01;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = m_t*U_sb_01*(m_t*IT_0000*IT_0010*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0002*IT_0059*U_sb_01*conj(U_st_01) +
       U_sb_11*(mu_h*conj(U_st_01) + m_t*IT_0010*IT_0059*conj(U_st_11)));
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0002*U_sb_11
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0005*((-0.5)*IT_0008*IT_0069 +
       m_W*IT_0000*IT_0002*U_sb_01*conj(U_st_01)) + A_t*IT_0000*U_sb_01*conj
      (U_st_11));
    const ccomplex_t IT_0071 = cpow((-2)*s_23 + IT_0020 + IT_0062 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0068*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0005*U_st_01;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = pow(m_d, 2);
    const ccomplex_t IT_0078 = A_t*IT_0000*conj(U_st_11) + (-1.4142135623731)
      *e_em*IT_0005*(m_W*IT_0000*IT_0002*conj(U_st_01) + (-0.5)*IT_0008*(IT_0002
      *IT_0059*IT_0077*conj(U_st_01) + m_t*(m_t*IT_0000*IT_0010*conj(U_st_01) +
       mu_h*conj(U_st_11))));
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*V_td*IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_23 + IT_0020 + IT_0062 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = IT_0076*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0005*conj(U_sb_00)*U_st_01;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = m_t*U_sb_00*(m_t*IT_0000*IT_0010*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0002*IT_0059*U_sb_00*conj(U_st_01) +
       U_sb_10*(mu_h*conj(U_st_01) + m_t*IT_0010*IT_0059*conj(U_st_11)));
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0002*U_sb_10
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0005*((-0.5)*IT_0008*IT_0086 +
       m_W*IT_0000*IT_0002*U_sb_00*conj(U_st_01)) + A_t*IT_0000*U_sb_00*conj
      (U_st_11));
    const ccomplex_t IT_0088 = cpow((-2)*s_23 + IT_0020 + IT_0062 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0085*IT_0090;
    const ccomplex_t IT_0092 = -IT_0066 + -IT_0074 + -IT_0083 + -IT_0091;
    const ccomplex_t IT_0093 = s_14*s_34*IT_0053;
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = s_13 + IT_0094;
    const ccomplex_t IT_0096 = -IT_0024 + -IT_0041 + -IT_0051 + IT_0066 +
       IT_0074 + IT_0083 + IT_0091;
    const ccomplex_t IT_0097 = s_14*s_24*IT_0053;
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = s_12 + IT_0098;
    const ccomplex_t IT_0100 = (-3)*conj(IT_0052);
    const ccomplex_t IT_0101 = pow(s_34, 2);
    const ccomplex_t IT_0102 = IT_0053*IT_0101;
    const ccomplex_t IT_0103 = s_24*s_34*IT_0053;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = s_23 + IT_0104;
    const ccomplex_t IT_0106 = pow(s_24, 2);
    const ccomplex_t IT_0107 = IT_0053*IT_0106;
    const ccomplex_t IT_0108 = (-3)*IT_0096;
    const ccomplex_t IT_0109 = (-3)*IT_0052*(conj(IT_0052)*(IT_0020 + -IT_0055
      ) + conj(IT_0092)*IT_0095 + conj(IT_0096)*IT_0099) + IT_0096*IT_0099
      *IT_0100 + (-3)*IT_0092*((-0.333333333333333)*IT_0095*IT_0100 + conj
      (IT_0092)*(IT_0062 + -IT_0102) + conj(IT_0096)*IT_0105) + (conj(IT_0092)
      *IT_0105 + conj(IT_0096)*(IT_0020 + -IT_0107))*IT_0108;
    return create_ccomplex_return(IT_0109);
}

