#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_tau_ss_R_to_tau_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_tau_ss_R_to_tau_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_tau, 2);
    const ccomplex_t IT_0011 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0012 = cpow(s_13 + (-0.5)*IT_0010 + (-0.5)*IT_0011 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em
      *m_tau*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = cpow(s_13 + (-0.5)*IT_0010 + (-0.5)*IT_0011 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0014 + (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0014 + IT_0022;
    const ccomplex_t IT_0025 = m_c*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = m_c*s_23;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0004;
    const ccomplex_t IT_0030 = IT_0018*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0004;
    const ccomplex_t IT_0033 = IT_0006*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0035 = IT_0012*IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (-0.5)*IT_0031 + (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = s_24*m_tau;
    const ccomplex_t IT_0038 = s_23*s_24;
    const ccomplex_t IT_0039 = pow(m_ss_R, 2);
    const ccomplex_t IT_0040 = s_34*IT_0039;
    const ccomplex_t IT_0041 = m_c*m_tau;
    const ccomplex_t IT_0042 = 6*s_34;
    const ccomplex_t IT_0043 = 6*conj(IT_0023)*(IT_0026*IT_0027 + -IT_0036
      *IT_0037 + 0.166666666666667*IT_0023*(12*IT_0038 + (-6)*IT_0040)) + 6*conj
      (IT_0026)*(IT_0023*IT_0027 + -IT_0036*IT_0041 + 0.166666666666667*IT_0026
      *IT_0042) + (-6)*conj(IT_0036)*(IT_0023*IT_0037 + IT_0026*IT_0041 + (
      -0.166666666666667)*IT_0036*IT_0042);
    return create_ccomplex_return(IT_0043);
}

