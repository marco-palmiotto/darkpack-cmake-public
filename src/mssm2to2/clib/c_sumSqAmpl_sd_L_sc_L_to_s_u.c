#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_sc_L_to_s_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_sc_L_to_s_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
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
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*conj(V_cs)*e_em*conj(V_Wp1)
      *IT_0001;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0006 = pow(m_s, 2);
    const ccomplex_t IT_0007 = pow(m_sc_L, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_23 + IT_0006 + IT_0007 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0004*IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*conj(V_cs)*e_em*conj(V_Wp2)
      *IT_0001;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0014 = cpow((-2)*s_23 + IT_0006 + IT_0007 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0013*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = sin(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0001*IT_0016*IT_0018;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0003*IT_0020;
    const ccomplex_t IT_0022 = IT_0008*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = m_u*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0001*IT_0016*IT_0018;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0011*IT_0026;
    const ccomplex_t IT_0028 = IT_0014*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = m_u*IT_0029;
    const ccomplex_t IT_0031 = -IT_0009 + -IT_0015 + IT_0024 + IT_0030;
    const ccomplex_t IT_0032 = -IT_0023 + -IT_0029;
    const ccomplex_t IT_0033 = m_u*s_13;
    const ccomplex_t IT_0034 = 18*s_34;
    const ccomplex_t IT_0035 = m_s*m_u;
    const ccomplex_t IT_0036 = cos(beta);
    const ccomplex_t IT_0037 = cpow(IT_0036, -1);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_cs)*e_em*IT_0001*IT_0016*IT_0037;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0020*IT_0039;
    const ccomplex_t IT_0041 = IT_0005*IT_0008*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_cs)*e_em*IT_0001*IT_0016*IT_0037;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = IT_0026*IT_0043;
    const ccomplex_t IT_0045 = IT_0013*IT_0014*IT_0044;
    const ccomplex_t IT_0046 = IT_0002*IT_0039;
    const ccomplex_t IT_0047 = IT_0008*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = m_u*IT_0048;
    const ccomplex_t IT_0050 = IT_0010*IT_0043;
    const ccomplex_t IT_0051 = IT_0014*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = m_u*IT_0052;
    const ccomplex_t IT_0054 = -IT_0041 + -IT_0045 + IT_0049 + IT_0053;
    const ccomplex_t IT_0055 = (-18)*conj(IT_0054);
    const ccomplex_t IT_0056 = m_s*s_14;
    const ccomplex_t IT_0057 = -IT_0048 + -IT_0052;
    const ccomplex_t IT_0058 = (-18)*conj(IT_0057);
    const ccomplex_t IT_0059 = pow(m_sd_L, 2);
    const ccomplex_t IT_0060 = s_34*IT_0059;
    const ccomplex_t IT_0061 = (-18)*IT_0060;
    const ccomplex_t IT_0062 = s_13*s_14;
    const ccomplex_t IT_0063 = 36*IT_0062;
    const ccomplex_t IT_0064 = IT_0061 + IT_0063;
    const ccomplex_t IT_0065 = m_s*m_u*IT_0059;
    const ccomplex_t IT_0066 = (-18)*conj(IT_0032);
    const ccomplex_t IT_0067 = IT_0031*(18*conj(IT_0032)*IT_0033 + conj
      (IT_0031)*IT_0034 + IT_0035*IT_0055 + IT_0056*IT_0058) + IT_0032*(18*conj
      (IT_0031)*IT_0033 + IT_0055*IT_0056 + conj(IT_0032)*IT_0064 + IT_0058
      *IT_0065) + IT_0054*((-18)*conj(IT_0031)*IT_0035 + IT_0034*conj(IT_0054) +
       18*IT_0033*conj(IT_0057) + IT_0056*IT_0066) + IT_0057*(18*IT_0033*conj
      (IT_0054) + (-18)*conj(IT_0031)*IT_0056 + conj(IT_0057)*IT_0064 + IT_0065
      *IT_0066);
    return create_ccomplex_return(IT_0067);
}

