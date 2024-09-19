#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sd_R_to_A0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sd_R_to_A0_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_d*e_em*mu_h*IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = A_b*V_td*conj(V_td);
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0006 + 0.5*IT_0010;
    const ccomplex_t IT_0012 = pow(m_A0, 2);
    const ccomplex_t IT_0013 = pow(m_sd_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0003*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0001;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*A_b*conj(V_td)*V_ts*IT_0007;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0026*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0007
      *U_sb_00;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = e_em*IT_0001;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0001*IT_0004;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = pow(m_su_L, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_12 + -IT_0013 + -IT_0039 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0036*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0026*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0007
      *U_sb_01;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0045*IT_0050;
    const ccomplex_t IT_0052 = -IT_0017 + -IT_0025 + -IT_0034 + IT_0043 + 
      -IT_0051;
    const ccomplex_t IT_0053 = pow(m_W, -2);
    const ccomplex_t IT_0054 = pow(s_34, 2);
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = IT_0017 + IT_0025 + IT_0034 + IT_0051;
    const ccomplex_t IT_0057 = s_24*s_34*IT_0053;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = s_23 + IT_0058;
    const ccomplex_t IT_0060 = pow(s_24, 2);
    const ccomplex_t IT_0061 = IT_0053*IT_0060;
    return create_ccomplex_return((-3)*IT_0052*(conj(IT_0052)*(IT_0012 + 
      -IT_0055) + conj(IT_0056)*IT_0059) + (-3)*IT_0056*(conj(IT_0052)*IT_0059 +
       conj(IT_0056)*(IT_0013 + -IT_0061)));
}

