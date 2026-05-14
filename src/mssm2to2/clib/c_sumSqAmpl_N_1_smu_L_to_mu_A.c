#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_smu_L_to_mu_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_smu_L_to_mu_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = conj(N_B1)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = conj(N_W1)*e_em;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0006 + IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_mu, 2);
    const ccomplex_t IT_0015 = pow(m_N_1, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0014 + IT_0015 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = pow(m_smu_L, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + IT_0014 + -IT_0015 + -IT_0020 
      + -reg_prop, -1);
    const ccomplex_t IT_0022 = -IT_0000;
    const ccomplex_t IT_0023 = IT_0013*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0019 + 2*IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = cos(beta);
    const ccomplex_t IT_0029 = cpow(IT_0028, -1);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0008*IT_0027*IT_0029;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0016*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0001*IT_0033;
    const ccomplex_t IT_0035 = IT_0022*IT_0031;
    const ccomplex_t IT_0036 = IT_0021*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -IT_0034 + (-2)*IT_0037;
    const ccomplex_t IT_0039 = (-2)*IT_0014;
    const ccomplex_t IT_0040 = IT_0000*IT_0018;
    const ccomplex_t IT_0041 = IT_0000*IT_0033;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = (-2)*s_13;
    const ccomplex_t IT_0044 = s_34*IT_0043;
    const ccomplex_t IT_0045 = -IT_0019;
    const ccomplex_t IT_0046 = -IT_0037;
    const ccomplex_t IT_0047 = m_mu*s_14*m_N_1;
    const ccomplex_t IT_0048 = s_13*s_14;
    const ccomplex_t IT_0049 = s_34*IT_0015;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0048 + IT_0050;
    const ccomplex_t IT_0052 = s_14*IT_0043;
    const ccomplex_t IT_0053 = IT_0015*IT_0043;
    const ccomplex_t IT_0054 = pow(s_13, 2);
    const ccomplex_t IT_0055 = (-2)*IT_0054;
    const ccomplex_t IT_0056 = m_mu*m_N_1;
    const ccomplex_t IT_0057 = (-2)*s_14;
    const ccomplex_t IT_0058 = (-2)*s_34;
    const ccomplex_t IT_0059 = (-2)*IT_0015;
    const ccomplex_t IT_0060 = (conj(IT_0026)*IT_0038 + IT_0026*conj(IT_0038))
      *IT_0039 + IT_0043*(conj(IT_0026)*IT_0034 + IT_0026*conj(IT_0034) + conj
      (IT_0038)*IT_0045 + IT_0038*conj(IT_0045)) + (conj(IT_0034)*IT_0040 +
       IT_0034*conj(IT_0040) + conj(IT_0042)*IT_0045 + IT_0042*conj(IT_0045))
      *IT_0057 + (conj(IT_0038)*IT_0040 + IT_0038*conj(IT_0040) + conj(IT_0026)
      *IT_0042 + IT_0026*conj(IT_0042))*IT_0058 + (conj(IT_0034)*IT_0045 +
       IT_0034*conj(IT_0045))*IT_0059;
    const ccomplex_t IT_0061 = m_mu*s_34*m_N_1;
    const ccomplex_t IT_0062 = (-2)*conj(IT_0026);
    const ccomplex_t IT_0063 = s_14*IT_0014;
    const ccomplex_t IT_0064 = (-2)*conj(IT_0038);
    const ccomplex_t IT_0065 = s_14*s_34;
    const ccomplex_t IT_0066 = (-2)*IT_0038;
    const ccomplex_t IT_0067 = (-0.5)*IT_0063;
    const ccomplex_t IT_0068 = (-2)*IT_0026;
    const ccomplex_t IT_0069 = (-0.5)*IT_0061;
    const ccomplex_t IT_0070 = s_13*(IT_0026*conj(IT_0026) + IT_0038*conj
      (IT_0038))*IT_0039 + (conj(IT_0026)*IT_0040 + IT_0026*conj(IT_0040) + conj
      (IT_0038)*IT_0042 + IT_0038*conj(IT_0042))*IT_0044 + (-2)*conj(IT_0045)*
      (IT_0046*IT_0047 + 2*IT_0025*IT_0051 + (-0.5)*IT_0040*IT_0052 + (-0.5)
      *IT_0045*IT_0053 + (-0.5)*IT_0026*IT_0055) + IT_0045*(conj(IT_0040)
      *IT_0052 + conj(IT_0026)*IT_0055) + (-2)*conj(IT_0034)*(IT_0025*IT_0047 +
       2*IT_0046*IT_0051 + (-0.5)*IT_0042*IT_0052 + (-0.5)*IT_0034*IT_0053 + (
      -0.5)*IT_0038*IT_0055) + IT_0034*(conj(IT_0042)*IT_0052 + conj(IT_0038)
      *IT_0055) + IT_0056*IT_0060 + IT_0025*(IT_0062*IT_0063 + IT_0061*IT_0064) 
      + IT_0046*(IT_0061*IT_0062 + IT_0063*IT_0064) + (-2)*conj(IT_0025)*
      (IT_0034*IT_0047 + 2*IT_0045*IT_0051 + (-4)*IT_0025*IT_0065 + IT_0067
      *IT_0068 + IT_0066*IT_0069) + (-2)*conj(IT_0046)*(IT_0045*IT_0047 + 2
      *IT_0034*IT_0051 + (-4)*IT_0046*IT_0065 + IT_0066*IT_0067 + IT_0068
      *IT_0069);
    return create_ccomplex_return(IT_0070);
}

