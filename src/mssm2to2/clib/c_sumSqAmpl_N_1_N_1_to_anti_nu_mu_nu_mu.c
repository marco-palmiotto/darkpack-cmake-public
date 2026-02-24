#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_anti_nu_mu_nu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_anti_nu_mu_nu_mu(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t IT_0000 = pow(m_N_1, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = e_em*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0004, -1);
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + IT_0009);
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = IT_0001*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0015 = IT_0005*IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0014;
    const ccomplex_t IT_0017 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0018 = IT_0005*IT_0017;
    const ccomplex_t IT_0019 = IT_0008*IT_0017;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0015 + IT_0016 + 
      -IT_0018 + -IT_0019);
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0020;
    const ccomplex_t IT_0025 = IT_0013*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = conj(N_B1)*e_em;
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = conj(N_W1)*e_em;
    const ccomplex_t IT_0031 = IT_0007*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0029 + -IT_0032);
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = N_B1*e_em;
    const ccomplex_t IT_0036 = IT_0003*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = N_W1*e_em;
    const ccomplex_t IT_0039 = IT_0007*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + -IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0034*IT_0042;
    const ccomplex_t IT_0044 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = s_34*IT_0000;
    const ccomplex_t IT_0049 = 4*IT_0048;
    const ccomplex_t IT_0050 = conj(IT_0047)*IT_0049;
    const ccomplex_t IT_0051 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0052 = IT_0043*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = s_13*s_24;
    const ccomplex_t IT_0056 = 8*IT_0055;
    const ccomplex_t IT_0057 = conj(IT_0054)*IT_0056;
    const ccomplex_t IT_0058 = conj(IT_0023) + conj(IT_0026);
    const ccomplex_t IT_0059 = 16*IT_0055;
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0022;
    const ccomplex_t IT_0062 = (-0.5)*IT_0025;
    const ccomplex_t IT_0063 = conj(IT_0061) + conj(IT_0062);
    const ccomplex_t IT_0064 = 8*IT_0048;
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = IT_0050 + IT_0057 + IT_0060 + IT_0065;
    const ccomplex_t IT_0067 = IT_0049*conj(IT_0054);
    const ccomplex_t IT_0068 = conj(IT_0023)*IT_0064;
    const ccomplex_t IT_0069 = s_14*s_23;
    const ccomplex_t IT_0070 = 16*IT_0069;
    const ccomplex_t IT_0071 = IT_0063*IT_0070;
    const ccomplex_t IT_0072 = conj(IT_0026)*IT_0064;
    const ccomplex_t IT_0073 = 8*IT_0069;
    const ccomplex_t IT_0074 = conj(IT_0047)*IT_0073;
    const ccomplex_t IT_0075 = IT_0067 + IT_0068 + IT_0071 + IT_0072 + IT_0074;
    const ccomplex_t IT_0076 = IT_0049*IT_0054;
    const ccomplex_t IT_0077 = 2*IT_0054*(conj(IT_0047)*IT_0048 + 2*conj
      (IT_0054)*IT_0055 + 0.5*IT_0056*IT_0058) + (IT_0023 + IT_0026)*IT_0066 + 4
      *IT_0047*(0.5*IT_0048*conj(IT_0054) + 0.25*IT_0049*IT_0058 + conj(IT_0047)
      *IT_0069 + 0.25*IT_0063*IT_0073) + (IT_0061 + IT_0062)*IT_0075 + IT_0063
      *IT_0076;
    return create_ccomplex_return(IT_0077);
}

