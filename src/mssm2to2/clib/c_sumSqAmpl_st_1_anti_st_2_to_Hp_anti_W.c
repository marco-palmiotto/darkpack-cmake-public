#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_Hp_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_Hp_anti_W(
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
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_st_1 = param->m_st_1;
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
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0001*conj(U_st_01);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = pow(m_d, 2);
    const ccomplex_t IT_0008 = cpow(IT_0004, -1);
    const ccomplex_t IT_0009 = cpow(IT_0005, -1);
    const ccomplex_t IT_0010 = A_t*IT_0004*U_st_10 + (-1.4142135623731)*e_em
      *IT_0001*(m_W*IT_0004*IT_0005*U_st_00 + (-0.5)*IT_0006*(IT_0005*IT_0007
      *IT_0008*U_st_00 + m_t*(m_t*IT_0004*IT_0009*U_st_00 + mu_h*U_st_10)));
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*conj(V_td)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_Hp, 2);
    const ccomplex_t IT_0013 = pow(m_st_1, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0003*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*U_sb_01*conj(U_st_01);
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = m_t*conj(U_sb_01)*(m_t*IT_0004*IT_0009*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0005*IT_0008*conj(U_sb_01)*U_st_00 + conj
      (U_sb_11)*(mu_h*U_st_00 + m_t*IT_0008*IT_0009*U_st_10));
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_11)*U_st_00 + (-1.4142135623731)*e_em*IT_0001*((-0.5)*IT_0006
      *IT_0020 + m_W*IT_0004*IT_0005*conj(U_sb_01)*U_st_00) + A_t*IT_0004*conj
      (U_sb_01)*U_st_10);
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0001*conj(U_st_01);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = pow(m_s, 2);
    const ccomplex_t IT_0029 = A_t*IT_0004*U_st_10 + (-1.4142135623731)*e_em
      *IT_0001*(m_W*IT_0004*IT_0005*U_st_00 + (-0.5)*IT_0006*(IT_0005*IT_0008
      *IT_0028*U_st_00 + m_t*(m_t*IT_0004*IT_0009*U_st_00 + mu_h*U_st_10)));
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*conj(V_ts)*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*U_sb_00*conj(U_st_01);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = m_t*conj(U_sb_00)*(m_t*IT_0004*IT_0009*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0005*IT_0008*conj(U_sb_00)*U_st_00 + conj
      (U_sb_10)*(mu_h*U_st_00 + m_t*IT_0008*IT_0009*U_st_10));
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_10)*U_st_00 + (-1.4142135623731)*e_em*IT_0001*((-0.5)*IT_0006
      *IT_0037 + m_W*IT_0004*IT_0005*conj(U_sb_00)*U_st_00) + A_t*IT_0004*conj
      (U_sb_00)*U_st_10);
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0036*IT_0041;
    const ccomplex_t IT_0043 = IT_0017 + IT_0025 + IT_0034 + IT_0042;
    const ccomplex_t IT_0044 = pow(m_st_2, 2);
    const ccomplex_t IT_0045 = pow(m_W, -2);
    const ccomplex_t IT_0046 = pow(s_24, 2);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = e_em*IT_0001;
    const ccomplex_t IT_0049 = m_t*e_em*mu_h*IT_0006*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0050 = IT_0001*IT_0049;
    const ccomplex_t IT_0051 = m_t*e_em*mu_h*IT_0006*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0052 = IT_0001*IT_0051;
    const ccomplex_t IT_0053 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0054 = IT_0004*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0057 = IT_0004*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (-0.5)*IT_0050 + 0.5*IT_0052 + (-0.5)*IT_0055 +
       0.5*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_12 + -IT_0013 + -IT_0044 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0048*IT_0062;
    const ccomplex_t IT_0064 = cos(alpha);
    const ccomplex_t IT_0065 = sin(alpha);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0005
      *IT_0064 + -IT_0004*IT_0065);
    const ccomplex_t IT_0067 = m_t*mu_h*IT_0006;
    const ccomplex_t IT_0068 = IT_0004*IT_0064;
    const ccomplex_t IT_0069 = -IT_0005*IT_0065 + IT_0068;
    const ccomplex_t IT_0070 = cos(theta_W);
    const ccomplex_t IT_0071 = cpow(IT_0070, -2);
    const ccomplex_t IT_0072 = pow(m_t, 2);
    const ccomplex_t IT_0073 = (-2)*U_st_00*((0.5*m_W*IT_0069*(IT_0001 + -1./3
      *IT_0000*IT_0071) + IT_0001*IT_0006*IT_0009*IT_0065*IT_0072)*conj(U_st_01)
       + -1./2*IT_0001*IT_0009*IT_0064*IT_0067*conj(U_st_11)) + U_st_10*(IT_0001
      *IT_0009*IT_0064*IT_0067*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0000
      *IT_0069*IT_0071 + 1.5*IT_0001*IT_0006*IT_0009*IT_0065*IT_0072)*conj
      (U_st_11));
    const ccomplex_t IT_0074 = (0 + _Complex_I*-1)*(e_em*IT_0073 + (
      -1.4142135623731)*A_t*IT_0065*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + -IT_0013 + -IT_0044 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0066*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0004
      *IT_0064 + IT_0005*IT_0065);
    const ccomplex_t IT_0081 = -IT_0080;
    const ccomplex_t IT_0082 = IT_0004*IT_0065;
    const ccomplex_t IT_0083 = IT_0005*IT_0064;
    const ccomplex_t IT_0084 = IT_0082 + IT_0083;
    const ccomplex_t IT_0085 = 2*U_st_00*((IT_0001*IT_0006*IT_0009*IT_0064
      *IT_0072 + (-0.5)*m_W*(IT_0001 + -1./3*IT_0000*IT_0071)*IT_0084)*conj
      (U_st_01) + 1./2*IT_0001*IT_0009*IT_0065*IT_0067*conj(U_st_11)) + U_st_10*
      (IT_0001*IT_0009*IT_0065*IT_0067*conj(U_st_01) + (-1.33333333333333)*((
      -1.5)*IT_0001*IT_0006*IT_0009*IT_0064*IT_0072 + m_W*IT_0000*IT_0071
      *IT_0084)*conj(U_st_11));
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(e_em*IT_0085 +
       1.4142135623731*A_t*IT_0064*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0013 + -IT_0044 + -reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0081*IT_0090;
    const ccomplex_t IT_0092 = IT_0063 + IT_0079 + IT_0091;
    const ccomplex_t IT_0093 = s_24*s_34*IT_0045;
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = s_23 + IT_0094;
    const ccomplex_t IT_0096 = pow(s_34, 2);
    const ccomplex_t IT_0097 = IT_0045*IT_0096;
    return create_ccomplex_return((-3)*IT_0043*(conj(IT_0043)*(IT_0044 + 
      -IT_0047) + conj(IT_0092)*IT_0095) + (-3)*IT_0092*(conj(IT_0043)*IT_0095 +
       conj(IT_0092)*(IT_0012 + -IT_0097)));
}

