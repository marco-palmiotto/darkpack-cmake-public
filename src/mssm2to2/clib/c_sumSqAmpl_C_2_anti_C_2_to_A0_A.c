#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_A0_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_A0_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = s_14*s_24;
    const ccomplex_t IT_0001 = pow(m_A0, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0001 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = V_u2*U_Wm2;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = U_d2*V_Wp2;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (IT_0006 + IT_0009)*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*e_em*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = IT_0003*IT_0014;
    const ccomplex_t IT_0016 = IT_0002*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = cpow((-2)*s_13 + IT_0001 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0025 = IT_0007*IT_0024;
    const ccomplex_t IT_0026 = IT_0011*(IT_0023 + IT_0025);
    const ccomplex_t IT_0027 = 1.4142135623731*e_em*IT_0026;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0003*IT_0028;
    const ccomplex_t IT_0030 = IT_0002*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0018*IT_0029;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = pow(m_C_2, 2);
    const ccomplex_t IT_0036 = s_14*IT_0035;
    const ccomplex_t IT_0037 = (-2)*IT_0031;
    const ccomplex_t IT_0038 = 2*conj(IT_0017);
    const ccomplex_t IT_0039 = (-2)*IT_0017;
    const ccomplex_t IT_0040 = 2*conj(IT_0031);
    const ccomplex_t IT_0041 = s_24*IT_0035;
    const ccomplex_t IT_0042 = 2*IT_0033;
    const ccomplex_t IT_0043 = 2*IT_0017;
    const ccomplex_t IT_0044 = 2*IT_0020;
    const ccomplex_t IT_0045 = 2*IT_0031;
    const ccomplex_t IT_0046 = (-2)*IT_0041*(conj(IT_0031)*IT_0037 + IT_0031
      *conj(IT_0037) + conj(IT_0017)*IT_0039 + IT_0017*conj(IT_0039) + (-0.5)
      *IT_0038*IT_0042 + (-0.5)*conj(IT_0042)*IT_0043 + (-0.5)*IT_0040*IT_0044 +
       (-0.5)*conj(IT_0044)*IT_0045);
    const ccomplex_t IT_0047 = s_12*s_24;
    const ccomplex_t IT_0048 = (-2)*IT_0047;
    const ccomplex_t IT_0049 = IT_0036 + IT_0048;
    const ccomplex_t IT_0050 = 2*s_12*IT_0035;
    const ccomplex_t IT_0051 = 2*cpow(IT_0035, 2);
    const ccomplex_t IT_0052 = s_12*s_14;
    const ccomplex_t IT_0053 = (-2)*IT_0052;
    const ccomplex_t IT_0054 = IT_0041 + IT_0053;
    const ccomplex_t IT_0055 = 2*IT_0054;
    const ccomplex_t IT_0056 = 2*IT_0036;
    const ccomplex_t IT_0057 = (-2)*s_12*IT_0035;
    const ccomplex_t IT_0058 = pow(s_12, 2);
    const ccomplex_t IT_0059 = (-2)*IT_0058;
    const ccomplex_t IT_0060 = 2*IT_0041;
    const ccomplex_t IT_0061 = 8*IT_0000;
    const ccomplex_t IT_0062 = (-2)*IT_0036;
    const ccomplex_t IT_0063 = (-8)*IT_0000*(conj(IT_0017)*IT_0021 + -IT_0017*
      (conj(IT_0017) + -conj(IT_0021)) + conj(IT_0031)*IT_0034 + -IT_0031*(conj
      (IT_0031) + -conj(IT_0034))) + IT_0036*(IT_0037*IT_0038 + IT_0039*IT_0040)
       + IT_0046 + (IT_0040*IT_0042 + IT_0038*IT_0044)*IT_0049 + conj(IT_0037)*
      (IT_0036*IT_0043 + IT_0044*IT_0050 + IT_0039*IT_0051 + IT_0034*IT_0055 +
       IT_0021*IT_0056 + IT_0037*IT_0057 + IT_0042*IT_0059) + conj(IT_0039)*
      (IT_0036*IT_0045 + IT_0042*IT_0050 + IT_0037*IT_0051 + IT_0021*IT_0055 +
       IT_0034*IT_0056 + IT_0039*IT_0057 + IT_0044*IT_0059) + conj(IT_0044)*
      (IT_0043*IT_0049 + IT_0037*IT_0050 + IT_0042*IT_0051 + IT_0044*IT_0057 +
       IT_0039*IT_0059 + IT_0034*IT_0060 + IT_0021*IT_0062) + conj(IT_0042)*
      (IT_0045*IT_0049 + IT_0039*IT_0050 + IT_0044*IT_0051 + IT_0042*IT_0057 +
       IT_0037*IT_0059 + IT_0021*IT_0060 + IT_0034*IT_0062) + conj(IT_0034)*
      (IT_0037*IT_0055 + IT_0039*IT_0056 + IT_0044*IT_0060 + IT_0034*IT_0061 +
       IT_0042*IT_0062) + conj(IT_0021)*(IT_0039*IT_0055 + IT_0037*IT_0056 +
       IT_0042*IT_0060 + IT_0021*IT_0061 + IT_0044*IT_0062);
    return create_ccomplex_return(IT_0063);
}

