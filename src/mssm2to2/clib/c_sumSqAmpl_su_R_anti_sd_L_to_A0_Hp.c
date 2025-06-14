#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_sd_L_to_A0_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_sd_L_to_A0_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t Gamma_ul = param->Gamma_ul;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_d*m_u*V_ud
      *e_em*IT_0000*IT_0002*IT_0004*IT_0006;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = m_d*e_em*mu_h*IT_0000;
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = A_b*V_td*conj(V_td);
    const ccomplex_t IT_0012 = IT_0003*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0010 + (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = pow(m_A0, 2);
    const ccomplex_t IT_0017 = pow(m_sd_L, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_s*m_u*V_us
      *e_em*IT_0000*IT_0002*IT_0004*IT_0006;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*A_b*conj(V_td)*V_ts*IT_0003;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = m_W*e_em*IT_0006;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *mu_h*IT_0000*IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = pow(m_W, 2);
    const ccomplex_t IT_0035 = pow(m_su_R, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + -IT_0017 + IT_0034 + -IT_0035 
      + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0033*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0000*IT_0006*IT_0039*(mu_h*conj(U_sb_00) + m_b*IT_0002*IT_0004*conj
      (U_sb_10))*V_ub_mod;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0003
      *U_sb_10;
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0000*IT_0006*IT_0039*(mu_h*conj(U_sb_01) + m_b*IT_0002*IT_0004*conj
      (U_sb_11))*V_ub_mod;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0003
      *U_sb_11;
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = m_u*e_em*mu_h*IT_0000*IT_0006;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = pow(m_d, 2);
    const ccomplex_t IT_0059 = pow(m_u, 2);
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0006*(IT_0000*IT_0002*IT_0003*IT_0058 + (-2)*IT_0001*(m_W*IT_0003 + (
      -0.5)*IT_0000*IT_0004*IT_0059));
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = IT_0057*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0016 + IT_0035 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0020 + IT_0028 + -IT_0038 + IT_0047 +
       IT_0055 + IT_0065;
    return create_ccomplex_return(3*IT_0066*conj(IT_0066));
}

