#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_snu_tau_to_anti_tau_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_snu_tau_to_anti_tau_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = s_34*m_C_1*m_tau;
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0002;
    const ccomplex_t IT_0004 = -IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0006 = -IT_0005;
    const ccomplex_t IT_0007 = IT_0004*IT_0006;
    const ccomplex_t IT_0008 = pow(m_C_1, 2);
    const ccomplex_t IT_0009 = pow(m_tau, 2);
    const ccomplex_t IT_0010 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0011 = cpow(s_12 + 0.5*IT_0008 + (-0.5)*IT_0009 + 0.5
      *IT_0010 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = s_14*IT_0009;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = cos(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0002*IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0006*IT_0020;
    const ccomplex_t IT_0022 = IT_0011*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0005*IT_0020;
    const ccomplex_t IT_0026 = cpow(s_23 + (-0.5)*IT_0009 + (-0.5)*IT_0010 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = s_13*s_34;
    const ccomplex_t IT_0031 = (-2)*IT_0030;
    const ccomplex_t IT_0032 = IT_0015 + IT_0031;
    const ccomplex_t IT_0033 = (-2)*conj(IT_0023);
    const ccomplex_t IT_0034 = s_34*IT_0008;
    const ccomplex_t IT_0035 = IT_0004*IT_0005;
    const ccomplex_t IT_0036 = IT_0026*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = s_14*m_C_1*m_tau;
    const ccomplex_t IT_0040 = (-2)*conj(IT_0037);
    const ccomplex_t IT_0041 = -IT_0013;
    const ccomplex_t IT_0042 = 2*IT_0000;
    const ccomplex_t IT_0043 = s_14*s_34;
    const ccomplex_t IT_0044 = 8*IT_0043;
    const ccomplex_t IT_0045 = IT_0024*IT_0044;
    const ccomplex_t IT_0046 = (-2)*IT_0023;
    const ccomplex_t IT_0047 = (-2)*IT_0037;
    const ccomplex_t IT_0048 = -IT_0028;
    const ccomplex_t IT_0049 = (-2)*IT_0034;
    const ccomplex_t IT_0050 = m_C_1*m_tau;
    const ccomplex_t IT_0051 = 2*s_13*IT_0050;
    const ccomplex_t IT_0052 = pow(s_13, 2);
    const ccomplex_t IT_0053 = 2*IT_0052;
    const ccomplex_t IT_0054 = 2*IT_0039;
    const ccomplex_t IT_0055 = s_13*s_14;
    const ccomplex_t IT_0056 = (-0.5)*IT_0034;
    const ccomplex_t IT_0057 = IT_0055 + IT_0056;
    const ccomplex_t IT_0058 = 4*IT_0057;
    const ccomplex_t IT_0059 = (-2)*s_13;
    const ccomplex_t IT_0060 = (-2)*IT_0008;
    const ccomplex_t IT_0061 = (-2)*IT_0039;
    const ccomplex_t IT_0062 = IT_0014 + IT_0038;
    const ccomplex_t IT_0063 = (-2)*IT_0015;
    const ccomplex_t IT_0064 = (-2)*IT_0009;
    const ccomplex_t IT_0065 = (-2)*IT_0000;
    const ccomplex_t IT_0066 = (-2)*IT_0032;
    const ccomplex_t IT_0067 = 0.5*IT_0041;
    const ccomplex_t IT_0068 = 0.5*IT_0048;
    const ccomplex_t IT_0069 = (IT_0000*IT_0014 + IT_0015*IT_0024 + IT_0029
      *IT_0032)*IT_0033 + (IT_0034*IT_0038 + IT_0029*IT_0039)*IT_0040 + conj
      (IT_0029)*(IT_0041*IT_0042 + IT_0029*IT_0044 + IT_0045 + IT_0032*IT_0046 +
       IT_0039*IT_0047 + IT_0048*IT_0049) + conj(IT_0048)*(IT_0029*IT_0049 +
       IT_0041*IT_0051 + IT_0023*IT_0053 + IT_0014*IT_0054 + IT_0024*IT_0058 +
       IT_0008*IT_0048*IT_0059 + IT_0037*IT_0050*IT_0060 + IT_0038*IT_0061) + 4
      *conj(IT_0014)*(0.25*IT_0000*IT_0046 + 0.25*IT_0048*IT_0054 + IT_0037
      *IT_0057 + 0.25*IT_0044*IT_0062 + 0.25*IT_0041*IT_0063) + conj(IT_0041)*
      (IT_0029*IT_0042 + IT_0048*IT_0051 + IT_0037*IT_0053 + IT_0009*IT_0041
      *IT_0059 + IT_0014*IT_0063 + IT_0023*IT_0050*IT_0064 + IT_0024*IT_0065 +
       IT_0038*IT_0066) + 2*conj(IT_0023)*(IT_0000*IT_0038 + 0.5*IT_0009*IT_0023
      *IT_0059 + IT_0050*(s_13*IT_0037 + IT_0064*IT_0067) + IT_0053*IT_0068) + 2
      *conj(IT_0024)*(IT_0037*IT_0039 + 0.5*IT_0029*IT_0044 + 0.5*IT_0045 + 0.5
      *IT_0015*IT_0046 + IT_0065*IT_0067 + IT_0058*IT_0068) + 2*conj(IT_0038)*
      (IT_0000*IT_0023 + 0.5*IT_0034*IT_0047 + 0.5*IT_0044*IT_0062 + IT_0066
      *IT_0067 + IT_0061*IT_0068) + 2*conj(IT_0037)*(IT_0024*IT_0039 + 2*IT_0014
      *IT_0057 + 0.5*IT_0008*IT_0037*IT_0059 + IT_0053*IT_0067 + IT_0050*(s_13
      *IT_0023 + IT_0060*IT_0068));
    return create_ccomplex_return(IT_0069);
}

