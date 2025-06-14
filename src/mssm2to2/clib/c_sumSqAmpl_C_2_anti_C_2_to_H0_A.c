#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_H0_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_H0_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = s_14*s_24;
    const ccomplex_t IT_0001 = pow(m_H0, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0001 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = U_d2*V_Wp2;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = V_u2*U_Wm2;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (IT_0006 + IT_0009)*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*e_em*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = IT_0003*IT_0015;
    const ccomplex_t IT_0017 = IT_0002*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0001 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0026 = IT_0007*IT_0025;
    const ccomplex_t IT_0027 = IT_0011*(IT_0024 + IT_0026);
    const ccomplex_t IT_0028 = 1.4142135623731*e_em*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = IT_0003*IT_0030;
    const ccomplex_t IT_0032 = IT_0002*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0019*IT_0031;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = pow(m_C_2, 2);
    const ccomplex_t IT_0038 = s_14*IT_0037;
    const ccomplex_t IT_0039 = (-2)*IT_0033;
    const ccomplex_t IT_0040 = 2*conj(IT_0018);
    const ccomplex_t IT_0041 = (-2)*IT_0018;
    const ccomplex_t IT_0042 = 2*conj(IT_0033);
    const ccomplex_t IT_0043 = s_24*IT_0037;
    const ccomplex_t IT_0044 = 2*IT_0035;
    const ccomplex_t IT_0045 = 2*IT_0018;
    const ccomplex_t IT_0046 = 2*IT_0021;
    const ccomplex_t IT_0047 = 2*IT_0033;
    const ccomplex_t IT_0048 = (-2)*IT_0043*(conj(IT_0033)*IT_0039 + IT_0033
      *conj(IT_0039) + conj(IT_0018)*IT_0041 + IT_0018*conj(IT_0041) + (-0.5)
      *IT_0040*IT_0044 + (-0.5)*conj(IT_0044)*IT_0045 + (-0.5)*IT_0042*IT_0046 +
       (-0.5)*conj(IT_0046)*IT_0047);
    const ccomplex_t IT_0049 = s_12*s_24;
    const ccomplex_t IT_0050 = (-2)*IT_0049;
    const ccomplex_t IT_0051 = IT_0038 + IT_0050;
    const ccomplex_t IT_0052 = 2*s_12*IT_0037;
    const ccomplex_t IT_0053 = 2*cpow(IT_0037, 2);
    const ccomplex_t IT_0054 = 2*IT_0038;
    const ccomplex_t IT_0055 = s_12*s_14;
    const ccomplex_t IT_0056 = (-0.5)*IT_0043;
    const ccomplex_t IT_0057 = IT_0055 + IT_0056;
    const ccomplex_t IT_0058 = (-4)*IT_0057;
    const ccomplex_t IT_0059 = (-2)*s_12*IT_0037;
    const ccomplex_t IT_0060 = pow(s_12, 2);
    const ccomplex_t IT_0061 = (-2)*IT_0060;
    const ccomplex_t IT_0062 = 2*IT_0043;
    const ccomplex_t IT_0063 = 8*IT_0000;
    const ccomplex_t IT_0064 = (-2)*IT_0038;
    const ccomplex_t IT_0065 = (-8)*IT_0000*(conj(IT_0018)*IT_0022 + -IT_0018*
      (conj(IT_0018) + -conj(IT_0022)) + conj(IT_0033)*IT_0036 + -IT_0033*(conj
      (IT_0033) + -conj(IT_0036))) + IT_0038*(IT_0039*IT_0040 + IT_0041*IT_0042)
       + IT_0048 + (IT_0042*IT_0044 + IT_0040*IT_0046)*IT_0051 + conj(IT_0039)*
      (IT_0038*IT_0045 + IT_0046*IT_0052 + IT_0041*IT_0053 + IT_0022*IT_0054 +
       IT_0036*IT_0058 + IT_0039*IT_0059 + IT_0044*IT_0061) + conj(IT_0041)*
      (IT_0038*IT_0047 + IT_0044*IT_0052 + IT_0039*IT_0053 + IT_0036*IT_0054 +
       IT_0022*IT_0058 + IT_0041*IT_0059 + IT_0046*IT_0061) + conj(IT_0046)*
      (IT_0045*IT_0051 + IT_0039*IT_0052 + IT_0044*IT_0053 + IT_0046*IT_0059 +
       IT_0041*IT_0061 + IT_0036*IT_0062 + IT_0022*IT_0064) + conj(IT_0044)*
      (IT_0047*IT_0051 + IT_0041*IT_0052 + IT_0046*IT_0053 + IT_0044*IT_0059 +
       IT_0039*IT_0061 + IT_0022*IT_0062 + IT_0036*IT_0064) + conj(IT_0036)*
      (IT_0041*IT_0054 + IT_0039*IT_0058 + IT_0046*IT_0062 + IT_0036*IT_0063 +
       IT_0044*IT_0064) + conj(IT_0022)*(IT_0039*IT_0054 + IT_0041*IT_0058 +
       IT_0044*IT_0062 + IT_0022*IT_0063 + IT_0046*IT_0064);
    return create_ccomplex_return(IT_0065);
}

