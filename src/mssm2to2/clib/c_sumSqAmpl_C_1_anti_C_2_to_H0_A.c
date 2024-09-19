#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_H0_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_H0_A(
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
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = U_d2*V_Wp1;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = sin(alpha);
    const ccomplex_t IT_0005 = V_u1*U_Wm2;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (IT_0003 + IT_0006)*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*e_em*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0000*IT_0012;
    const ccomplex_t IT_0014 = pow(m_H0, 2);
    const ccomplex_t IT_0015 = pow(m_C_1, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0014 + IT_0015 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = 2*IT_0018;
    const ccomplex_t IT_0020 = -IT_0018;
    const ccomplex_t IT_0021 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = IT_0008*(IT_0022 + IT_0024);
    const ccomplex_t IT_0026 = 1.4142135623731*e_em*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0000*IT_0028;
    const ccomplex_t IT_0030 = IT_0016*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = 2*IT_0031;
    const ccomplex_t IT_0033 = -IT_0031;
    const ccomplex_t IT_0034 = pow(m_C_2, 2);
    const ccomplex_t IT_0035 = s_14*IT_0034;
    const ccomplex_t IT_0036 = s_14*s_24;
    const ccomplex_t IT_0037 = cpow((-2)*s_23 + IT_0014 + IT_0034 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0013*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = IT_0029*IT_0037;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = -conj(IT_0041);
    const ccomplex_t IT_0043 = conj(IT_0033) + IT_0042;
    const ccomplex_t IT_0044 = m_C_1*m_C_2;
    const ccomplex_t IT_0045 = (-2)*IT_0039;
    const ccomplex_t IT_0046 = (-2)*IT_0041;
    const ccomplex_t IT_0047 = 2*s_12;
    const ccomplex_t IT_0048 = 2*IT_0034;
    const ccomplex_t IT_0049 = IT_0015*(conj(IT_0045)*IT_0046 + IT_0045*conj
      (IT_0046)) + 0.5*(conj(IT_0032)*IT_0045 + IT_0032*conj(IT_0045) + conj
      (IT_0019)*IT_0046 + IT_0019*conj(IT_0046))*IT_0047 + 0.5*(conj(IT_0019)
      *IT_0032 + IT_0019*conj(IT_0032))*IT_0048;
    const ccomplex_t IT_0050 = s_14*m_C_1*m_C_2;
    const ccomplex_t IT_0051 = 2*IT_0045;
    const ccomplex_t IT_0052 = IT_0033 + IT_0041;
    const ccomplex_t IT_0053 = 2*conj(IT_0045);
    const ccomplex_t IT_0054 = 2*IT_0046;
    const ccomplex_t IT_0055 = 2*conj(IT_0046);
    const ccomplex_t IT_0056 = s_24*m_C_1*m_C_2;
    const ccomplex_t IT_0057 = 2*IT_0019;
    const ccomplex_t IT_0058 = 2*conj(IT_0019);
    const ccomplex_t IT_0059 = 2*IT_0032;
    const ccomplex_t IT_0060 = 2*conj(IT_0032);
    const ccomplex_t IT_0061 = s_12*s_24;
    const ccomplex_t IT_0062 = (-2)*IT_0061;
    const ccomplex_t IT_0063 = IT_0035 + IT_0062;
    const ccomplex_t IT_0064 = (-2)*s_12;
    const ccomplex_t IT_0065 = IT_0034*IT_0064;
    const ccomplex_t IT_0066 = s_24*IT_0015;
    const ccomplex_t IT_0067 = s_12*s_14;
    const ccomplex_t IT_0068 = (-0.5)*IT_0066;
    const ccomplex_t IT_0069 = IT_0067 + IT_0068;
    const ccomplex_t IT_0070 = pow(s_12, 2);
    const ccomplex_t IT_0071 = (-2)*IT_0070;
    const ccomplex_t IT_0072 = (-2)*(conj(IT_0019)*IT_0020 + IT_0019*conj
      (IT_0020) + conj(IT_0032)*IT_0033 + IT_0032*conj(IT_0033))*IT_0035 + 8
      *IT_0036*((IT_0020 + -IT_0039)*(conj(IT_0020) + -conj(IT_0039)) + (IT_0033
       + -IT_0041)*IT_0043) + 2*IT_0044*IT_0049 + IT_0050*((conj(IT_0033) + conj
      (IT_0041))*IT_0051 + IT_0052*IT_0053 + (conj(IT_0020) + conj(IT_0039))
      *IT_0054 + (IT_0020 + IT_0039)*IT_0055) + IT_0056*((conj(IT_0033) + conj
      (IT_0041))*IT_0057 + IT_0052*IT_0058 + (conj(IT_0020) + conj(IT_0039))
      *IT_0059 + (IT_0020 + IT_0039)*IT_0060) + (conj(IT_0039)*IT_0057 + IT_0039
      *IT_0058 + conj(IT_0041)*IT_0059 + IT_0041*IT_0060)*IT_0063 + (IT_0019
      *conj(IT_0019) + IT_0032*conj(IT_0032))*IT_0065 + (-2)*conj(IT_0045)*((
      -0.5)*IT_0015*IT_0045*IT_0064 + IT_0039*IT_0066 + 2*IT_0020*IT_0069 + (
      -0.5)*IT_0019*IT_0071) + (-2)*IT_0045*(conj(IT_0039)*IT_0066 + 2*conj
      (IT_0020)*IT_0069 + (-0.5)*conj(IT_0019)*IT_0071) + (-2)*conj(IT_0046)*((
      -0.5)*IT_0015*IT_0046*IT_0064 + IT_0041*IT_0066 + 2*IT_0033*IT_0069 + (
      -0.5)*IT_0032*IT_0071) + (-2)*IT_0046*(conj(IT_0041)*IT_0066 + 2*conj
      (IT_0033)*IT_0069 + (-0.5)*conj(IT_0032)*IT_0071);
    return create_ccomplex_return(IT_0072);
}

