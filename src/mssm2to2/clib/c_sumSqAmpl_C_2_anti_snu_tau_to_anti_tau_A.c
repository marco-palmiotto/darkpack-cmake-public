#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_snu_tau_to_anti_tau_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_snu_tau_to_anti_tau_A(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_tau, 2);
    const ccomplex_t IT_0001 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0002 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = pow(m_C_2, 2);
    const ccomplex_t IT_0016 = cpow(s_12 + (-0.5)*IT_0000 + 0.5*IT_0001 + 0.5
      *IT_0015 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0017 = -IT_0003;
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = s_13*s_14;
    const ccomplex_t IT_0023 = s_34*IT_0015;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0022 + IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0008;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = IT_0002*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0017*IT_0027;
    const ccomplex_t IT_0033 = IT_0016*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = -IT_0034;
    const ccomplex_t IT_0037 = pow(s_13, 2);
    const ccomplex_t IT_0038 = 2*IT_0037;
    const ccomplex_t IT_0039 = s_14*s_34;
    const ccomplex_t IT_0040 = 0.5*IT_0013;
    const ccomplex_t IT_0041 = 8*IT_0021;
    const ccomplex_t IT_0042 = 8*IT_0035;
    const ccomplex_t IT_0043 = (-2)*s_13;
    const ccomplex_t IT_0044 = IT_0000*IT_0043;
    const ccomplex_t IT_0045 = m_C_2*m_tau;
    const ccomplex_t IT_0046 = 2*s_13;
    const ccomplex_t IT_0047 = (-2)*IT_0015;
    const ccomplex_t IT_0048 = (-2)*IT_0000;
    const ccomplex_t IT_0049 = (-2)*IT_0014;
    const ccomplex_t IT_0050 = (-2)*conj(IT_0014);
    const ccomplex_t IT_0051 = s_14*m_C_2*m_tau;
    const ccomplex_t IT_0052 = (-2)*IT_0040;
    const ccomplex_t IT_0053 = (-2)*conj(IT_0040);
    const ccomplex_t IT_0054 = s_14*IT_0000;
    const ccomplex_t IT_0055 = s_13*s_34;
    const ccomplex_t IT_0056 = (-2)*IT_0055;
    const ccomplex_t IT_0057 = IT_0054 + IT_0056;
    const ccomplex_t IT_0058 = (-2)*IT_0036;
    const ccomplex_t IT_0059 = (-2)*conj(IT_0036);
    const ccomplex_t IT_0060 = s_34*m_C_2*m_tau;
    const ccomplex_t IT_0061 = (-0.5)*IT_0057;
    const ccomplex_t IT_0062 = 4*IT_0014*conj(IT_0021)*IT_0025 + (-2)*conj
      (IT_0030)*(IT_0023*IT_0031 + (-2)*IT_0025*IT_0035 + (-0.5)*IT_0036*IT_0038
      ) + 8*IT_0039*((conj(IT_0021) + conj(IT_0040))*(IT_0040 + 0.125*IT_0041) +
       (conj(IT_0031) + conj(IT_0035))*(IT_0031 + 0.125*IT_0042)) + 4*conj
      (IT_0014)*(IT_0021*IT_0025 + 0.25*IT_0014*IT_0015*IT_0043) + (-2)*IT_0030*
      (IT_0023*conj(IT_0031) + (-2)*IT_0025*conj(IT_0035) + (-0.5)*conj(IT_0036)
      *IT_0038 + (-0.5)*IT_0015*conj(IT_0030)*IT_0043) + IT_0036*conj(IT_0036)
      *IT_0044 + IT_0045*((conj(IT_0020)*IT_0030 + IT_0020*conj(IT_0030) + conj
      (IT_0014)*IT_0036 + IT_0014*conj(IT_0036))*IT_0046 + (conj(IT_0014)
      *IT_0030 + IT_0014*conj(IT_0030))*IT_0047 + (conj(IT_0020)*IT_0036 +
       IT_0020*conj(IT_0036))*IT_0048) + IT_0023*(conj(IT_0040)*IT_0049 +
       IT_0040*IT_0050) + 2*IT_0051*(conj(IT_0014)*IT_0035 + IT_0014*conj
      (IT_0035) + 0.5*conj(IT_0031)*IT_0049 + 0.5*IT_0031*IT_0050 + conj(IT_0030
      )*(IT_0021 + 0.5*IT_0052) + IT_0030*(conj(IT_0021) + 0.5*IT_0053)) +
       IT_0057*(conj(IT_0031)*IT_0058 + IT_0031*IT_0059) + IT_0054*(conj(IT_0035
      )*IT_0058 + IT_0035*IT_0059) + 2*(conj(IT_0020)*(IT_0031 + -IT_0035) +
       IT_0020*(conj(IT_0031) + -conj(IT_0035)) + conj(IT_0036)*IT_0040 +
       IT_0036*conj(IT_0040) + 0.5*conj(IT_0021)*IT_0058 + 0.5*IT_0021*IT_0059)
      *IT_0060 + (-2)*conj(IT_0020)*((-0.5)*IT_0014*IT_0038 + IT_0021*IT_0054 +
       IT_0052*IT_0061) + (-2)*IT_0020*((-0.5)*conj(IT_0014)*IT_0038 + (-0.5)
      *conj(IT_0020)*IT_0044 + conj(IT_0021)*IT_0054 + IT_0053*IT_0061);
    return create_ccomplex_return(IT_0062);
}

