#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sd_L_to_A0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sd_L_to_A0_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
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
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0000
      *IT_0002*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0004 = -IT_0003;
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0005
      *U_sb_10;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = pow(m_A0, 2);
    const ccomplex_t IT_0009 = pow(m_sd_L, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = e_em*IT_0002;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = pow(m_d, 2);
    const ccomplex_t IT_0018 = cos(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = pow(m_u, 2);
    const ccomplex_t IT_0021 = cpow(IT_0005, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0002*(IT_0005*IT_0016*IT_0017*IT_0019 + (-2)*IT_0018*(m_W*IT_0005 + (
      -0.5)*IT_0016*IT_0020*IT_0021));
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = pow(m_su_L, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + -IT_0009 + -IT_0024 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0015*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0000
      *IT_0002*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0005
      *U_sb_11;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = -IT_0013 + IT_0028 + -IT_0036;
    const ccomplex_t IT_0038 = pow(m_W, -2);
    const ccomplex_t IT_0039 = pow(s_34, 2);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = IT_0013 + IT_0036;
    const ccomplex_t IT_0042 = s_24*s_34*IT_0038;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = s_23 + IT_0043;
    const ccomplex_t IT_0045 = pow(s_24, 2);
    const ccomplex_t IT_0046 = IT_0038*IT_0045;
    return create_ccomplex_return((-3)*IT_0037*(conj(IT_0037)*(IT_0008 + 
      -IT_0040) + conj(IT_0041)*IT_0044) + (-3)*IT_0041*(conj(IT_0037)*IT_0044 +
       conj(IT_0041)*(IT_0009 + -IT_0046)));
}

