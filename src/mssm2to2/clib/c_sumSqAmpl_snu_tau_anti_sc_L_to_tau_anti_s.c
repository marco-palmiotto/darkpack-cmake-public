#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_tau_anti_sc_L_to_tau_anti_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_tau_anti_sc_L_to_tau_anti_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
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
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0001;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = pow(m_tau, 2);
    const ccomplex_t IT_0007 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0008 = cpow(s_13 + (-0.5)*IT_0006 + (-0.5)*IT_0007 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp2*IT_0001;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0001;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = cpow(s_13 + (-0.5)*IT_0006 + (-0.5)*IT_0007 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = (-0.5)*IT_0010 + (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = cos(beta);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0001*IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0012*IT_0024;
    const ccomplex_t IT_0026 = IT_0015*IT_0019*IT_0025;
    const ccomplex_t IT_0027 = m_tau*IT_0010;
    const ccomplex_t IT_0028 = m_tau*IT_0017;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em
      *m_tau*IT_0001*IT_0020*IT_0022;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = IT_0003*IT_0031;
    const ccomplex_t IT_0033 = IT_0008*IT_0029*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0026 + 0.5*IT_0027 + 0.5*IT_0028 + (
      -0.5)*IT_0033;
    const ccomplex_t IT_0035 = s_14*m_tau;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0001*IT_0020*IT_0022;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0004*IT_0037;
    const ccomplex_t IT_0039 = IT_0008*IT_0029*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0001*IT_0020*IT_0022;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = IT_0013*IT_0041;
    const ccomplex_t IT_0043 = IT_0015*IT_0019*IT_0042;
    const ccomplex_t IT_0044 = IT_0031*IT_0037;
    const ccomplex_t IT_0045 = IT_0008*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = m_tau*IT_0046;
    const ccomplex_t IT_0048 = IT_0024*IT_0041;
    const ccomplex_t IT_0049 = IT_0015*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = m_tau*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0039 + (-0.5)*IT_0043 + 0.5*IT_0047 +
       0.5*IT_0051;
    const ccomplex_t IT_0053 = m_s*s_13;
    const ccomplex_t IT_0054 = (-0.5)*IT_0046 + (-0.5)*IT_0050;
    const ccomplex_t IT_0055 = m_s*m_tau*IT_0007;
    const ccomplex_t IT_0056 = s_34*IT_0007;
    const ccomplex_t IT_0057 = (-6)*IT_0056;
    const ccomplex_t IT_0058 = s_13*s_14;
    const ccomplex_t IT_0059 = 12*IT_0058;
    const ccomplex_t IT_0060 = IT_0057 + IT_0059;
    const ccomplex_t IT_0061 = 6*s_34;
    const ccomplex_t IT_0062 = m_s*m_tau;
    const ccomplex_t IT_0063 = (-6)*conj(IT_0034);
    const ccomplex_t IT_0064 = (-6)*IT_0034;
    const ccomplex_t IT_0065 = -IT_0018;
    const ccomplex_t IT_0066 = 6*conj(IT_0018)*(IT_0034*IT_0035 + -IT_0052
      *IT_0053 + -IT_0054*IT_0055 + 0.166666666666667*IT_0018*IT_0060) + 6*conj
      (IT_0034)*(IT_0018*IT_0035 + 0.166666666666667*IT_0034*IT_0061) + (IT_0053
      *IT_0054 + IT_0052*IT_0062)*IT_0063 + 6*conj(IT_0052)*(IT_0035*IT_0054 +
       0.166666666666667*IT_0052*IT_0061 + 0.166666666666667*IT_0062*IT_0064 +
       IT_0053*IT_0065) + 6*conj(IT_0054)*(IT_0035*IT_0052 + 0.166666666666667
      *IT_0054*IT_0060 + 0.166666666666667*IT_0053*IT_0064 + IT_0055*IT_0065);
    return create_ccomplex_return(IT_0066);
}

