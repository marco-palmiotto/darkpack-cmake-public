#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_L_to_A0_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_L_to_A0_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = pow(m_W, -2);
    const ccomplex_t IT_0003 = pow(m_s, 2);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -2);
    const ccomplex_t IT_0006 = pow(m_c, 2);
    const ccomplex_t IT_0007 = cpow(IT_0001, -2);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = 1.4142135623731*conj(V_cs)*IT_0000*(IT_0001*
      (IT_0001 + -IT_0001*IT_0002*IT_0003*IT_0005) + -IT_0004*(IT_0004 + 
      -IT_0002*IT_0004*IT_0006*IT_0007))*IT_0009;
    const ccomplex_t IT_0011 = (-0.25)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = cpow(IT_0008, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = m_s*e_em*mu_h*IT_0012;
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = A_b*V_ts*conj(V_ts);
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0017 + (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0015*IT_0021;
    const ccomplex_t IT_0023 = pow(m_A0, 2);
    const ccomplex_t IT_0024 = pow(m_ss_L, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0023 + IT_0024 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = m_W*e_em*IT_0013;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = cpow(IT_0004, 2);
    const ccomplex_t IT_0031 = cpow(IT_0001, 2);
    const ccomplex_t IT_0032 = (0 + _Complex_I*(-1.4142135623731))*conj(V_cs)
      *e_em*IT_0013*((IT_0003 + -IT_0006)*IT_0012 + -m_W*(IT_0030 + -IT_0031));
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0029*IT_0033;
    const ccomplex_t IT_0035 = pow(m_W, 2);
    const ccomplex_t IT_0036 = pow(m_sc_L, 2);
    const ccomplex_t IT_0037 = cpow((-2)*s_12 + -IT_0024 + IT_0035 + -IT_0036 
      + -reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0034*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0001, -1);
    const ccomplex_t IT_0041 = cpow(IT_0004, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0013*(m_W*IT_0001*IT_0004*conj(U_sb_01) + (-0.5)*IT_0012*(IT_0004
      *IT_0006*IT_0040*conj(U_sb_01) + m_b*(m_b*IT_0001*IT_0041*conj(U_sb_01) +
       mu_h*conj(U_sb_11))));
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0001
      *U_sb_11;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0023 + IT_0024 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0024 + -IT_0036 + -reg_prop, -1);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0013;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = e_em*IT_0013;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0024 + 0.5*IT_0036;
    const ccomplex_t IT_0057 = s_12 + IT_0056;
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = IT_0052*IT_0053;
    const ccomplex_t IT_0060 = (-0.5)*s_13 + (-0.5)*s_23 + 0.5*IT_0024 + 0.5
      *IT_0036;
    const ccomplex_t IT_0061 = s_12 + IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = -IT_0051;
    const ccomplex_t IT_0064 = IT_0053*IT_0063;
    const ccomplex_t IT_0065 = s_12 + IT_0024;
    const ccomplex_t IT_0066 = -s_23;
    const ccomplex_t IT_0067 = IT_0065 + IT_0066;
    const ccomplex_t IT_0068 = IT_0064*IT_0067;
    const ccomplex_t IT_0069 = IT_0054*IT_0063;
    const ccomplex_t IT_0070 = IT_0065*IT_0069;
    const ccomplex_t IT_0071 = -IT_0062 + 0.5*IT_0068 + (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0058 + IT_0071;
    const ccomplex_t IT_0073 = IT_0050*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = m_c*e_em*mu_h*IT_0012*IT_0013;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*mu_h*IT_0012*IT_0013;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_13 + IT_0023 + IT_0036 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_cd)
      *e_em*mu_h*IT_0012*IT_0013;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*A_b*V_td*conj(V_ts)*IT_0001;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_23 + IT_0023 + IT_0024 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_23 + IT_0023 + IT_0024 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0092 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0001
      *U_sb_10;
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0013*(m_W*IT_0001*IT_0004*conj(U_sb_00) + (-0.5)*IT_0012*(IT_0004
      *IT_0006*IT_0040*conj(U_sb_00) + m_b*(m_b*IT_0001*IT_0041*conj(U_sb_00) +
       mu_h*conj(U_sb_10))));
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0093*IT_0095;
    const ccomplex_t IT_0097 = IT_0091*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0011 + IT_0027 + -IT_0039 + IT_0049 + (-2)
      *IT_0074 + IT_0082 + IT_0090 + IT_0098;
    return create_ccomplex_return(3*IT_0099*conj(IT_0099));
}

