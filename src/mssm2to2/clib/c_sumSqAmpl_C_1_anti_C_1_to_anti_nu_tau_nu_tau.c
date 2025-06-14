#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_anti_nu_tau_nu_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_anti_nu_tau_nu_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = e_em*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0002, -1);
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(IT_0004 + IT_0007);
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0006*IT_0010;
    const ccomplex_t IT_0013 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + -IT_0012 + (-2)
      *IT_0014);
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_1, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = s_13*s_24;
    const ccomplex_t IT_0024 = e_em*U_Wm1*conj(U_stau_01);
    const ccomplex_t IT_0025 = IT_0005*IT_0024;
    const ccomplex_t IT_0026 = cos(beta);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = IT_0005*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = U_d1*e_em*m_tau*IT_0029*conj(U_stau_11);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0025 + (-0.5)*IT_0032);
    const ccomplex_t IT_0034 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0035 = IT_0005*IT_0034;
    const ccomplex_t IT_0036 = conj(U_d1)*e_em*m_tau*IT_0029*U_stau_11;
    const ccomplex_t IT_0037 = IT_0028*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + (-0.5)*IT_0038);
    const ccomplex_t IT_0040 = IT_0033*IT_0039;
    const ccomplex_t IT_0041 = cpow(s_23 + (-0.5)*IT_0018 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = e_em*U_Wm1*conj(U_stau_00);
    const ccomplex_t IT_0045 = IT_0005*IT_0044;
    const ccomplex_t IT_0046 = U_d1*e_em*m_tau*IT_0029*conj(U_stau_10);
    const ccomplex_t IT_0047 = IT_0028*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + (-0.5)*IT_0048);
    const ccomplex_t IT_0050 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0051 = IT_0005*IT_0050;
    const ccomplex_t IT_0052 = conj(U_d1)*e_em*m_tau*IT_0029*U_stau_10;
    const ccomplex_t IT_0053 = IT_0028*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0051 + (-0.5)*IT_0054);
    const ccomplex_t IT_0056 = IT_0049*IT_0055;
    const ccomplex_t IT_0057 = cpow(s_23 + (-0.5)*IT_0018 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = (-0.5)*IT_0043 + (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0062 = IT_0003*IT_0061;
    const ccomplex_t IT_0063 = IT_0006*IT_0061;
    const ccomplex_t IT_0064 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0065 = IT_0006*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0062 + -IT_0063 + (-2)
      *IT_0065);
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0009*IT_0067;
    const ccomplex_t IT_0069 = IT_0019*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = s_14*s_23;
    const ccomplex_t IT_0072 = conj(IT_0060) + 2*conj(IT_0070);
    const ccomplex_t IT_0073 = s_34*IT_0018;
    return create_ccomplex_return(16*IT_0022*conj(IT_0022)*IT_0023 + 4*
      (IT_0060 + 2*IT_0070)*IT_0071*IT_0072 + 4*(conj(IT_0022)*(IT_0060 + 2
      *IT_0070) + IT_0022*IT_0072)*IT_0073);
}

