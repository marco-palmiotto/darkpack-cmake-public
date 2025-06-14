#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_L_to_A0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_L_to_A0_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0001*conj(U_sb_00);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0004
      *U_sb_10;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_A0, 2);
    const ccomplex_t IT_0008 = pow(m_ss_L, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0003*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0001*conj(U_sb_01);
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0004
      *U_sb_11;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0014*IT_0019;
    const ccomplex_t IT_0021 = IT_0012 + IT_0020;
    const ccomplex_t IT_0022 = pow(m_W, -2);
    const ccomplex_t IT_0023 = pow(s_24, 2);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = e_em*IT_0001;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = pow(m_s, 2);
    const ccomplex_t IT_0029 = cos(beta);
    const ccomplex_t IT_0030 = cpow(IT_0029, -1);
    const ccomplex_t IT_0031 = pow(m_c, 2);
    const ccomplex_t IT_0032 = cpow(IT_0004, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0001*(IT_0004*IT_0027*IT_0028*IT_0030 + (-2)*IT_0029*(m_W*IT_0004
       + (-0.5)*IT_0027*IT_0031*IT_0032));
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = pow(m_sc_L, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + -IT_0008 + -IT_0035 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0026*IT_0038;
    const ccomplex_t IT_0040 = -IT_0012 + -IT_0020 + IT_0039;
    const ccomplex_t IT_0041 = s_24*s_34*IT_0022;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = s_23 + IT_0042;
    const ccomplex_t IT_0044 = pow(s_34, 2);
    const ccomplex_t IT_0045 = IT_0022*IT_0044;
    return create_ccomplex_return((-3)*IT_0021*(conj(IT_0021)*(IT_0008 + 
      -IT_0024) + conj(IT_0040)*IT_0043) + (-3)*IT_0040*(conj(IT_0021)*IT_0043 +
       conj(IT_0040)*(IT_0007 + -IT_0045)));
}

