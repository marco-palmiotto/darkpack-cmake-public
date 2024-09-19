#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_sb_2_to_s_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_sb_2_to_s_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_sc_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u2)*e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0003*IT_0005*IT_0007*U_sb_01;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = IT_0002*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u1)*e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0003*IT_0005*IT_0007*U_sb_01;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = IT_0015*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0014 + IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0025 = cos(beta);
    const ccomplex_t IT_0026 = cpow(IT_0025, -1);
    const ccomplex_t IT_0027 = IT_0007*IT_0026;
    const ccomplex_t IT_0028 = m_b*conj(U_d1)*V_tb*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0032 = IT_0007*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0030 + (-2)*IT_0032);
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = IT_0017*IT_0034;
    const ccomplex_t IT_0036 = IT_0015*IT_0024*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0038 = m_b*conj(U_d2)*V_tb*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0039 = IT_0027*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0042 = IT_0007*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0040 + (-2)*IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0009*IT_0044;
    const ccomplex_t IT_0046 = IT_0002*IT_0037*IT_0045;
    const ccomplex_t IT_0047 = IT_0036 + IT_0046;
    const ccomplex_t IT_0048 = m_t*s_13;
    const ccomplex_t IT_0049 = s_13*s_14;
    const ccomplex_t IT_0050 = s_34*IT_0001;
    const ccomplex_t IT_0051 = m_s*s_14;
    const ccomplex_t IT_0052 = m_s*IT_0023;
    const ccomplex_t IT_0053 = m_s*m_t;
    const ccomplex_t IT_0054 = 18*s_34;
    const ccomplex_t IT_0055 = (-18)*IT_0023;
    const ccomplex_t IT_0056 = conj(IT_0023)*((-18)*IT_0047*IT_0048 + IT_0023*
      (36*IT_0049 + (-18)*IT_0050) + (-18)*IT_0051*IT_0052) + conj(IT_0047)*(18
      *IT_0052*IT_0053 + IT_0047*IT_0054 + IT_0048*IT_0055) + conj(IT_0052)*(18
      *IT_0047*IT_0053 + IT_0052*IT_0054 + IT_0051*IT_0055);
    return create_ccomplex_return(IT_0056);
}

