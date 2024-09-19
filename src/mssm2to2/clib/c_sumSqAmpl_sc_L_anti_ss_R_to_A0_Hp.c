#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_R_to_A0_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_R_to_A0_Hp(
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
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_sl = param->Gamma_sl;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = pow(m_d, 2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = pow(m_c, 2);
    const ccomplex_t IT_0006 = cpow(IT_0004, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*conj(V_cd)
      *e_em*(IT_0000*IT_0001*IT_0003*IT_0004 + (-2)*IT_0002*(m_W*IT_0004 + (-0.5
      )*IT_0000*IT_0005*IT_0006))*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*A_b*V_td*conj(V_ts)*IT_0004;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = pow(m_A0, 2);
    const ccomplex_t IT_0015 = pow(m_ss_R, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_23 + IT_0014 + IT_0015 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = pow(m_s, 2);
    const ccomplex_t IT_0020 = (0 + _Complex_I*(-2.82842712474619))*conj(V_cs)
      *e_em*IT_0008*(IT_0002*(m_W*IT_0004 + (-0.5)*IT_0000*IT_0005*IT_0006) + 
      -1./2*IT_0000*IT_0003*IT_0004*IT_0019);
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = m_s*e_em*mu_h*IT_0000;
    const ccomplex_t IT_0023 = IT_0008*IT_0022;
    const ccomplex_t IT_0024 = A_b*V_ts*conj(V_ts);
    const ccomplex_t IT_0025 = IT_0004*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = 0.5*IT_0023 + 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0021*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0014 + IT_0015 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = m_W*e_em*IT_0008;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0000*IT_0003*IT_0004*IT_0008;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = pow(m_W, 2);
    const ccomplex_t IT_0038 = pow(m_sc_L, 2);
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + -IT_0015 + IT_0037 + -IT_0038 
      + -reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0036*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0008*(m_W*IT_0002*IT_0004*conj(U_sb_01) + (-0.5)*IT_0000*(IT_0002
      *IT_0005*IT_0006*conj(U_sb_01) + m_b*(m_b*IT_0003*IT_0004*conj(U_sb_01) +
       mu_h*conj(U_sb_11))));
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0004
      *U_sb_01;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0014 + IT_0015 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = m_c*e_em*mu_h*IT_0000*IT_0008;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1.4142135623731)*m_c*m_s*conj
      (V_cs)*e_em*IT_0000*IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_13 + IT_0014 + IT_0038 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0014 + IT_0015 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0059 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0004
      *U_sb_00;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0008*(m_W*IT_0002*IT_0004*conj(U_sb_00) + (-0.5)*IT_0000*(IT_0002
      *IT_0005*IT_0006*conj(U_sb_00) + m_b*(m_b*IT_0003*IT_0004*conj(U_sb_00) +
       mu_h*conj(U_sb_10))));
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = IT_0058*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0018 + IT_0031 + -IT_0041 + IT_0049 +
       IT_0057 + IT_0065;
    return create_ccomplex_return(3*IT_0066*conj(IT_0066));
}

