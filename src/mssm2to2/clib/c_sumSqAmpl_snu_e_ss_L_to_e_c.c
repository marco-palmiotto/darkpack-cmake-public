#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_e_ss_L_to_e_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_e_ss_L_to_e_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_e*U_d1*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = pow(m_e, 2);
    const ccomplex_t IT_0010 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_e*U_d2*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0013 + IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = cpow(IT_0023, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0001*IT_0003*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0007*IT_0026;
    const ccomplex_t IT_0028 = IT_0011*IT_0022*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0001*IT_0003*IT_0024;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = IT_0016*IT_0031;
    const ccomplex_t IT_0033 = IT_0018*IT_0029*IT_0032;
    const ccomplex_t IT_0034 = m_c*IT_0013;
    const ccomplex_t IT_0035 = m_c*IT_0020;
    const ccomplex_t IT_0036 = IT_0028 + IT_0033 + -IT_0034 + -IT_0035;
    const ccomplex_t IT_0037 = m_c*s_23;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0001;
    const ccomplex_t IT_0039 = IT_0002*IT_0038;
    const ccomplex_t IT_0040 = IT_0011*IT_0022*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0001;
    const ccomplex_t IT_0042 = IT_0014*IT_0041;
    const ccomplex_t IT_0043 = IT_0018*IT_0029*IT_0042;
    const ccomplex_t IT_0044 = IT_0026*IT_0038;
    const ccomplex_t IT_0045 = IT_0011*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = m_c*IT_0046;
    const ccomplex_t IT_0048 = IT_0031*IT_0041;
    const ccomplex_t IT_0049 = IT_0018*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = m_c*IT_0050;
    const ccomplex_t IT_0052 = IT_0040 + IT_0043 + -IT_0047 + -IT_0051;
    const ccomplex_t IT_0053 = m_e*s_24;
    const ccomplex_t IT_0054 = IT_0046 + IT_0050;
    const ccomplex_t IT_0055 = pow(m_ss_L, 2);
    const ccomplex_t IT_0056 = m_c*m_e*IT_0055;
    const ccomplex_t IT_0057 = s_34*IT_0055;
    const ccomplex_t IT_0058 = (-6)*IT_0057;
    const ccomplex_t IT_0059 = s_23*s_24;
    const ccomplex_t IT_0060 = 12*IT_0059;
    const ccomplex_t IT_0061 = IT_0058 + IT_0060;
    const ccomplex_t IT_0062 = m_c*m_e;
    const ccomplex_t IT_0063 = 6*s_34;
    const ccomplex_t IT_0064 = -conj(IT_0021);
    const ccomplex_t IT_0065 = (-6)*conj(IT_0036);
    const ccomplex_t IT_0066 = 6*IT_0021*(conj(IT_0036)*IT_0037 + -conj
      (IT_0052)*IT_0053 + -conj(IT_0054)*IT_0056 + 0.166666666666667*conj
      (IT_0021)*IT_0061) + (-6)*IT_0036*(IT_0053*conj(IT_0054) + conj(IT_0052)
      *IT_0062 + (-0.166666666666667)*conj(IT_0036)*IT_0063 + IT_0037*IT_0064) +
       6*IT_0054*(IT_0037*conj(IT_0052) + 0.166666666666667*conj(IT_0054)
      *IT_0061 + IT_0056*IT_0064 + 0.166666666666667*IT_0053*IT_0065) + 6
      *IT_0052*(IT_0037*conj(IT_0054) + 0.166666666666667*conj(IT_0052)*IT_0063 
      + IT_0053*IT_0064 + 0.166666666666667*IT_0062*IT_0065);
    return create_ccomplex_return(IT_0066);
}

