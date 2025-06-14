#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_sb_1_to_nu_tau_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_sb_1_to_nu_tau_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u1)*e_em*IT_0000*IT_0002*IT_0004*conj(U_sb_00);
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0004*IT_0010;
    const ccomplex_t IT_0012 = conj(U_d1)*e_em*m_tau*IT_0000*U_stau_10;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0008 + (-0.5)*IT_0014);
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0018 = pow(m_stau_1, 2);
    const ccomplex_t IT_0019 = cpow(s_13 + (-0.5)*IT_0018 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0016*IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0000*IT_0002*IT_0004*conj(U_sb_00);
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = conj(U_d2)*e_em*m_tau*IT_0000*U_stau_10;
    const ccomplex_t IT_0026 = IT_0011*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + (-0.5)*IT_0027);
    const ccomplex_t IT_0029 = IT_0022*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0031 = cpow(s_13 + (-0.5)*IT_0018 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0029*IT_0030*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0020 + 0.5*IT_0032;
    const ccomplex_t IT_0034 = V_tb*e_em*U_Wm1*conj(U_sb_00);
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = m_b*U_d1*V_tb*e_em*IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0037 = IT_0011*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + (-0.5)*IT_0038);
    const ccomplex_t IT_0040 = IT_0015*IT_0039;
    const ccomplex_t IT_0041 = IT_0019*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = V_tb*e_em*U_Wm2*conj(U_sb_00);
    const ccomplex_t IT_0044 = IT_0004*IT_0043;
    const ccomplex_t IT_0045 = m_b*U_d2*V_tb*e_em*IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0046 = IT_0011*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0044 + (-0.5)*IT_0047);
    const ccomplex_t IT_0049 = IT_0028*IT_0048;
    const ccomplex_t IT_0050 = IT_0031*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = 0.5*IT_0042 + 0.5*IT_0051;
    const ccomplex_t IT_0053 = m_t*s_13;
    const ccomplex_t IT_0054 = s_13*s_14;
    const ccomplex_t IT_0055 = s_34*IT_0018;
    return create_ccomplex_return(6*IT_0033*(s_34*conj(IT_0033) + -conj
      (IT_0052)*IT_0053) + (-6)*IT_0052*(conj(IT_0033)*IT_0053 + (
      -0.166666666666667)*conj(IT_0052)*(12*IT_0054 + (-6)*IT_0055)));
}

