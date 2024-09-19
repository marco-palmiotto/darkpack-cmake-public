#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_A0_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_A0_A(
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
    const creal_t m_C_1 = param->m_C_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = pow(m_A0, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = V_u1*U_Wm1;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = U_d1*V_Wp1;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (IT_0005 + IT_0008)*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*e_em*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0002*IT_0013;
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0022 = IT_0003*IT_0021;
    const ccomplex_t IT_0023 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0024 = IT_0006*IT_0023;
    const ccomplex_t IT_0025 = IT_0010*(IT_0022 + IT_0024);
    const ccomplex_t IT_0026 = 1.4142135623731*e_em*IT_0025;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = IT_0002*IT_0027;
    const ccomplex_t IT_0029 = IT_0001*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = IT_0018*IT_0028;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = s_14*s_24;
    const ccomplex_t IT_0035 = pow(m_C_1, 2);
    const ccomplex_t IT_0036 = 2*IT_0016;
    const ccomplex_t IT_0037 = (-2)*IT_0033;
    const ccomplex_t IT_0038 = (-2)*IT_0020;
    const ccomplex_t IT_0039 = 2*IT_0030;
    const ccomplex_t IT_0040 = 2*s_12;
    const ccomplex_t IT_0041 = 2*IT_0035;
    const ccomplex_t IT_0042 = s_24*IT_0035;
    const ccomplex_t IT_0043 = s_12*s_14;
    const ccomplex_t IT_0044 = (-2)*IT_0043;
    const ccomplex_t IT_0045 = IT_0042 + IT_0044;
    const ccomplex_t IT_0046 = 2*IT_0038;
    const ccomplex_t IT_0047 = 2*conj(IT_0038);
    const ccomplex_t IT_0048 = 2*IT_0037;
    const ccomplex_t IT_0049 = 2*conj(IT_0037);
    const ccomplex_t IT_0050 = s_14*IT_0035;
    const ccomplex_t IT_0051 = conj(IT_0033)*IT_0046 + (-2)*conj(IT_0031)*
      (IT_0039 + (-0.5)*IT_0046) + IT_0033*IT_0047 + (-2)*IT_0031*(conj(IT_0039)
       + (-0.5)*IT_0047) + conj(IT_0020)*IT_0048 + (-2)*conj(IT_0017)*(IT_0036 +
       (-0.5)*IT_0048) + IT_0020*IT_0049 + (-2)*IT_0017*(conj(IT_0036) + (-0.5)
      *IT_0049);
    const ccomplex_t IT_0052 = s_12*s_24;
    const ccomplex_t IT_0053 = (-2)*IT_0052;
    const ccomplex_t IT_0054 = IT_0050 + IT_0053;
    const ccomplex_t IT_0055 = 2*IT_0036;
    const ccomplex_t IT_0056 = 2*conj(IT_0036);
    const ccomplex_t IT_0057 = 2*IT_0039;
    const ccomplex_t IT_0058 = 2*conj(IT_0039);
    const ccomplex_t IT_0059 = conj(IT_0031)*IT_0055 + (-2)*conj(IT_0033)*
      (IT_0037 + (-0.5)*IT_0055) + IT_0031*IT_0056 + (-2)*IT_0033*(conj(IT_0037)
       + (-0.5)*IT_0056) + conj(IT_0017)*IT_0057 + (-2)*conj(IT_0020)*(IT_0038 +
       (-0.5)*IT_0057) + IT_0017*IT_0058 + (-2)*IT_0020*(conj(IT_0038) + (-0.5)
      *IT_0058);
    const ccomplex_t IT_0060 = (-2)*s_12*IT_0035;
    const ccomplex_t IT_0061 = pow(s_12, 2);
    const ccomplex_t IT_0062 = (-2)*IT_0061;
    const ccomplex_t IT_0063 = 8*((IT_0017 + -IT_0020)*(conj(IT_0017) + -conj
      (IT_0020)) + (IT_0031 + -IT_0033)*(conj(IT_0031) + -conj(IT_0033)))
      *IT_0034 + IT_0035*((conj(IT_0036)*IT_0037 + IT_0036*conj(IT_0037) + conj
      (IT_0038)*IT_0039 + IT_0038*conj(IT_0039))*IT_0040 + (conj(IT_0037)
      *IT_0038 + IT_0037*conj(IT_0038) + conj(IT_0036)*IT_0039 + IT_0036*conj
      (IT_0039))*IT_0041) + IT_0045*(conj(IT_0017)*IT_0046 + IT_0017*IT_0047 +
       conj(IT_0031)*IT_0048 + IT_0031*IT_0049) + IT_0050*IT_0051 + IT_0054*
      (conj(IT_0020)*IT_0055 + IT_0020*IT_0056 + conj(IT_0033)*IT_0057 + IT_0033
      *IT_0058) + IT_0042*IT_0059 + (IT_0036*conj(IT_0036) + IT_0037*conj
      (IT_0037) + IT_0038*conj(IT_0038) + IT_0039*conj(IT_0039))*IT_0060 + (conj
      (IT_0036)*IT_0038 + IT_0036*conj(IT_0038) + conj(IT_0037)*IT_0039 +
       IT_0037*conj(IT_0039))*IT_0062;
    return create_ccomplex_return(IT_0063);
}

