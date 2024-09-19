#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_snu_mu_to_anti_mu_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_snu_mu_to_anti_mu_A(
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
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0002;
    const ccomplex_t IT_0004 = -IT_0003;
    const ccomplex_t IT_0005 = IT_0000*IT_0004;
    const ccomplex_t IT_0006 = pow(m_mu, 2);
    const ccomplex_t IT_0007 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_23 + IT_0006 + IT_0007 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = pow(m_C_1, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_12 + IT_0006 + -IT_0007 + -IT_0011 
      + -reg_prop, -1);
    const ccomplex_t IT_0013 = -IT_0000;
    const ccomplex_t IT_0014 = IT_0004*IT_0013;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = s_14*s_34;
    const ccomplex_t IT_0018 = m_mu*s_14*m_C_1;
    const ccomplex_t IT_0019 = pow(m_W, -1);
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0002*IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0000*IT_0023;
    const ccomplex_t IT_0025 = IT_0008*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = 2*IT_0026;
    const ccomplex_t IT_0028 = m_mu*s_34*m_C_1;
    const ccomplex_t IT_0029 = IT_0013*IT_0023;
    const ccomplex_t IT_0030 = IT_0012*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = (-2)*IT_0031;
    const ccomplex_t IT_0033 = s_34*IT_0011;
    const ccomplex_t IT_0034 = (-2)*IT_0010;
    const ccomplex_t IT_0035 = 2*IT_0016;
    const ccomplex_t IT_0036 = s_14*IT_0006;
    const ccomplex_t IT_0037 = s_13*s_34;
    const ccomplex_t IT_0038 = (-2)*IT_0037;
    const ccomplex_t IT_0039 = IT_0036 + IT_0038;
    const ccomplex_t IT_0040 = s_13*s_14;
    const ccomplex_t IT_0041 = (-2)*IT_0040;
    const ccomplex_t IT_0042 = IT_0033 + IT_0041;
    const ccomplex_t IT_0043 = (-2)*IT_0016;
    const ccomplex_t IT_0044 = m_mu*m_C_1;
    const ccomplex_t IT_0045 = 2*s_13*IT_0044;
    const ccomplex_t IT_0046 = pow(s_13, 2);
    const ccomplex_t IT_0047 = 2*IT_0046;
    const ccomplex_t IT_0048 = -IT_0031;
    const ccomplex_t IT_0049 = 2*IT_0018;
    const ccomplex_t IT_0050 = (-2)*s_13*IT_0011;
    const ccomplex_t IT_0051 = (-2)*IT_0011*IT_0044;
    const ccomplex_t IT_0052 = (-2)*conj(IT_0010);
    const ccomplex_t IT_0053 = (-2)*conj(IT_0016);
    const ccomplex_t IT_0054 = -IT_0026;
    const ccomplex_t IT_0055 = (-2)*IT_0018;
    const ccomplex_t IT_0056 = 2*IT_0028;
    const ccomplex_t IT_0057 = (-2)*IT_0006;
    const ccomplex_t IT_0058 = IT_0044*IT_0057;
    const ccomplex_t IT_0059 = (-2)*IT_0028;
    const ccomplex_t IT_0060 = 0.5*conj(IT_0035);
    const ccomplex_t IT_0061 = (-2)*IT_0042;
    const ccomplex_t IT_0062 = 0.5*conj(IT_0048);
    const ccomplex_t IT_0063 = 0.5*conj(IT_0027);
    const ccomplex_t IT_0064 = (-2)*IT_0033;
    const ccomplex_t IT_0065 = 0.5*conj(IT_0054);
    const ccomplex_t IT_0066 = conj(IT_0032)*IT_0047 + conj(IT_0034)*IT_0051 +
       2*IT_0018*(conj(IT_0016) + 0.5*IT_0052) + 2*IT_0045*IT_0060 + 2*IT_0061
      *IT_0062 + 2*IT_0050*IT_0063 + 2*IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (-2)*IT_0036;
    const ccomplex_t IT_0068 = (-2)*IT_0039;
    const ccomplex_t IT_0069 = conj(IT_0034)*IT_0045 + 2*IT_0028*(conj(IT_0010
      ) + 0.5*IT_0053) + s_13*conj(IT_0032)*IT_0057 + 2*IT_0058*IT_0060 + 2
      *IT_0047*IT_0063 + 2*IT_0062*IT_0067 + 2*IT_0065*IT_0068;
    const ccomplex_t IT_0070 = conj(IT_0048) + conj(IT_0054);
    const ccomplex_t IT_0071 = 8*IT_0017;
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = 8*(IT_0010 + IT_0016)*(conj(IT_0010) + conj
      (IT_0016))*IT_0017 + 2*IT_0016*IT_0018*conj(IT_0027) + 2*IT_0010*IT_0028
      *conj(IT_0032) + IT_0034*(IT_0018*conj(IT_0027) + IT_0033*conj(IT_0034) +
       conj(IT_0035)*IT_0039) + (IT_0028*conj(IT_0032) + conj(IT_0035)*IT_0036 +
       conj(IT_0034)*IT_0042)*IT_0043 + IT_0034*(conj(IT_0032)*IT_0045 + conj
      (IT_0035)*IT_0047 + conj(IT_0048)*IT_0049 + conj(IT_0034)*IT_0050 + conj
      (IT_0027)*IT_0051 + IT_0033*IT_0052 + IT_0042*IT_0053 + conj(IT_0054)
      *IT_0055) + IT_0035*(conj(IT_0027)*IT_0045 + conj(IT_0034)*IT_0047 +
       IT_0039*IT_0052 + IT_0036*IT_0053 + conj(IT_0054)*IT_0056 + s_13*conj
      (IT_0035)*IT_0057 + conj(IT_0032)*IT_0058 + conj(IT_0048)*IT_0059) +
       IT_0027*IT_0066 + IT_0032*IT_0069 + IT_0048*(conj(IT_0034)*IT_0049 + conj
      (IT_0035)*IT_0059 + conj(IT_0027)*IT_0061 + conj(IT_0032)*IT_0067 +
       IT_0072) + IT_0054*(conj(IT_0034)*IT_0055 + conj(IT_0035)*IT_0056 + conj
      (IT_0027)*IT_0064 + conj(IT_0032)*IT_0068 + IT_0072);
    return create_ccomplex_return(IT_0073);
}

