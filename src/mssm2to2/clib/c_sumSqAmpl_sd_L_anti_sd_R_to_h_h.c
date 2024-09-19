#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_sd_R_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_sd_R_to_h_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = sin(alpha);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0002 = 0.5*IT_0001;
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0000*U_sb_01;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = IT_0002*IT_0004;
    const ccomplex_t IT_0006 = pow(m_h, 2);
    const ccomplex_t IT_0007 = pow(m_sd_L, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, -1);
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = cos(alpha);
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0000 + m_d*e_em*mu_h*IT_0011*IT_0013*IT_0014*IT_0016);
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = pow(m_d, 2);
    const ccomplex_t IT_0020 = sin(beta);
    const ccomplex_t IT_0021 = cos(theta_W);
    const ccomplex_t IT_0022 = cpow(IT_0021, -2);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0011
      *IT_0013*IT_0016*IT_0019 + (-0.5)*m_W*(IT_0000*IT_0012 + IT_0014*IT_0020)*
      (IT_0016 + 1./3*IT_0015*IT_0022));
    const ccomplex_t IT_0024 = IT_0018*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0006 + (0 + _Complex_I*1)
      *m_sd_L*Gamma_dl + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = pow(m_sd_R, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0006 + IT_0028 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0024*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*-3)*e_em*(IT_0000*IT_0011
      *IT_0013*IT_0016*IT_0019 + -1./3*m_W*IT_0015*(IT_0000*IT_0012 + IT_0014
      *IT_0020)*IT_0022);
    const ccomplex_t IT_0033 = (-0.333333333333333)*IT_0032;
    const ccomplex_t IT_0034 = IT_0018*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0006 + (0 + _Complex_I*1)
      *m_sd_R*Gamma_dr + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0034*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0000*U_sb_00;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = IT_0041*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0006 + IT_0028 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0046*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0006 + IT_0028 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0005*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0014 + -m_d*e_em*mu_h*IT_0000*IT_0011*IT_0013*IT_0016);
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0015*IT_0022;
    const ccomplex_t IT_0058 = IT_0016 + IT_0057;
    const ccomplex_t IT_0059 = cpow(IT_0000, 3);
    const ccomplex_t IT_0060 = cpow(IT_0014, 2);
    const ccomplex_t IT_0061 = cpow(IT_0000, 2);
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_W*e_em*IT_0058*(IT_0020
      *IT_0059 + (-5)*IT_0014*(IT_0000*IT_0014*IT_0020 + (-0.2)*IT_0012*(IT_0060
       + (-5)*IT_0061)));
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = IT_0056*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_12 + -IT_0007 + -IT_0028 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0007 + -IT_0028 + -reg_prop, -1);
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*m_W*e_em*IT_0058*(IT_0012
      *IT_0059 + -IT_0014*(IT_0014*(IT_0000*IT_0012 + IT_0014*IT_0020) + 
      -IT_0020*IT_0061));
    const ccomplex_t IT_0070 = 1.5*IT_0069;
    const ccomplex_t IT_0071 = IT_0018*IT_0070;
    const ccomplex_t IT_0072 = IT_0068*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0010 + IT_0027 + IT_0031 + IT_0037 + IT_0040
       + IT_0048 + IT_0051 + IT_0054 + -IT_0067 + -IT_0073;
    return create_ccomplex_return(3*IT_0074*conj(IT_0074));
}

