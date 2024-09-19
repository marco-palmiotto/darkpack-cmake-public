#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_su_L_to_nu_tau_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_su_L_to_nu_tau_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp1)*IT_0001;
    const ccomplex_t IT_0003 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = conj(U_d1)*e_em*m_tau*IT_0008*U_stau_10;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0004 + (-0.5)*IT_0011);
    const ccomplex_t IT_0013 = IT_0002*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0015 = pow(m_stau_1, 2);
    const ccomplex_t IT_0016 = cpow(s_13 + (-0.5)*IT_0015 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0013*IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp2)*IT_0001;
    const ccomplex_t IT_0019 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0020 = IT_0001*IT_0019;
    const ccomplex_t IT_0021 = conj(U_d2)*e_em*m_tau*IT_0008*U_stau_10;
    const ccomplex_t IT_0022 = IT_0007*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0020 + (-0.5)*IT_0023);
    const ccomplex_t IT_0025 = IT_0018*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0027 = cpow(s_13 + (-0.5)*IT_0015 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0025*IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (-0.5)*IT_0017 + (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*V_ud
      *e_em*IT_0001*IT_0006*IT_0008;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = IT_0012*IT_0031;
    const ccomplex_t IT_0033 = IT_0016*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*V_ud
      *e_em*IT_0001*IT_0006*IT_0008;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = IT_0024*IT_0036;
    const ccomplex_t IT_0038 = IT_0027*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = IT_0034 + IT_0039;
    const ccomplex_t IT_0041 = m_d*IT_0040;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = (-0.5)*IT_0034 + (-0.5)*IT_0039;
    const ccomplex_t IT_0044 = m_d*s_23;
    const ccomplex_t IT_0045 = s_23*s_24;
    const ccomplex_t IT_0046 = pow(m_su_L, 2);
    const ccomplex_t IT_0047 = s_34*IT_0046;
    const ccomplex_t IT_0048 = (-6)*s_34;
    const ccomplex_t IT_0049 = 6*s_34*(IT_0029*conj(IT_0029) + IT_0042*conj
      (IT_0042)) + 6*((conj(IT_0029) + -conj(IT_0042))*IT_0043 + (IT_0029 + 
      -IT_0042)*conj(IT_0043))*IT_0044 + IT_0043*conj(IT_0043)*(12*IT_0045 + (-6
      )*IT_0047) + (conj(IT_0029)*IT_0042 + IT_0029*conj(IT_0042))*IT_0048;
    return create_ccomplex_return(IT_0049);
}

