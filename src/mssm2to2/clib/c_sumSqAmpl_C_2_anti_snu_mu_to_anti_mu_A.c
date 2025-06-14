#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_snu_mu_to_anti_mu_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_snu_mu_to_anti_mu_A(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_mu, 2);
    const ccomplex_t IT_0001 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_mu*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = 2*IT_0013;
    const ccomplex_t IT_0015 = pow(m_C_2, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0015 
      + -reg_prop, -1);
    const ccomplex_t IT_0017 = -IT_0003;
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (-2)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0008;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = IT_0003*IT_0023;
    const ccomplex_t IT_0025 = IT_0002*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = (-2)*IT_0026;
    const ccomplex_t IT_0028 = IT_0017*IT_0023;
    const ccomplex_t IT_0029 = IT_0016*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = 2*IT_0030;
    const ccomplex_t IT_0032 = pow(s_13, 2);
    const ccomplex_t IT_0033 = 2*IT_0032;
    const ccomplex_t IT_0034 = s_14*s_34;
    const ccomplex_t IT_0035 = -IT_0013;
    const ccomplex_t IT_0036 = -IT_0020;
    const ccomplex_t IT_0037 = 8*IT_0026;
    const ccomplex_t IT_0038 = (-2)*IT_0000;
    const ccomplex_t IT_0039 = (-2)*s_13*IT_0015;
    const ccomplex_t IT_0040 = m_mu*m_C_2;
    const ccomplex_t IT_0041 = 2*s_13;
    const ccomplex_t IT_0042 = (-2)*IT_0015;
    const ccomplex_t IT_0043 = s_34*IT_0015;
    const ccomplex_t IT_0044 = (-2)*IT_0027;
    const ccomplex_t IT_0045 = (-2)*conj(IT_0027);
    const ccomplex_t IT_0046 = (-2)*IT_0014;
    const ccomplex_t IT_0047 = (-2)*conj(IT_0014);
    const ccomplex_t IT_0048 = s_13*s_14;
    const ccomplex_t IT_0049 = (-2)*IT_0048;
    const ccomplex_t IT_0050 = IT_0043 + IT_0049;
    const ccomplex_t IT_0051 = m_mu*s_14*m_C_2;
    const ccomplex_t IT_0052 = 2*(conj(IT_0014)*IT_0030 + IT_0014*conj(IT_0030
      ) + conj(IT_0027)*IT_0036 + IT_0027*conj(IT_0036) + 0.5*conj(IT_0035)
      *IT_0044 + 0.5*IT_0035*IT_0045 + 0.5*conj(IT_0026)*IT_0046 + 0.5*IT_0026
      *IT_0047)*IT_0051;
    const ccomplex_t IT_0053 = m_mu*s_34*m_C_2;
    const ccomplex_t IT_0054 = (-2)*IT_0021;
    const ccomplex_t IT_0055 = (-2)*conj(IT_0021);
    const ccomplex_t IT_0056 = (-2)*IT_0031;
    const ccomplex_t IT_0057 = (-2)*conj(IT_0031);
    const ccomplex_t IT_0058 = 2*IT_0053*(conj(IT_0021)*IT_0026 + IT_0021*conj
      (IT_0026) + conj(IT_0031)*IT_0035 + IT_0031*conj(IT_0035) + 0.5*conj
      (IT_0030)*IT_0054 + 0.5*IT_0030*IT_0055 + 0.5*conj(IT_0036)*IT_0056 + 0.5
      *IT_0036*IT_0057);
    const ccomplex_t IT_0059 = s_14*IT_0000;
    const ccomplex_t IT_0060 = s_13*s_34;
    const ccomplex_t IT_0061 = (-2)*IT_0060;
    const ccomplex_t IT_0062 = IT_0059 + IT_0061;
    const ccomplex_t IT_0063 = (conj(IT_0014)*IT_0021 + IT_0014*conj(IT_0021) 
      + conj(IT_0027)*IT_0031 + IT_0027*conj(IT_0031))*IT_0033 + 8*IT_0034*(
      (IT_0035 + IT_0036)*(conj(IT_0035) + conj(IT_0036)) + (conj(IT_0026) +
       conj(IT_0030))*(IT_0030 + 0.125*IT_0037)) + s_13*(IT_0021*conj(IT_0021) +
       IT_0031*conj(IT_0031))*IT_0038 + (IT_0014*conj(IT_0014) + IT_0027*conj
      (IT_0027))*IT_0039 + IT_0040*((conj(IT_0021)*IT_0031 + IT_0021*conj
      (IT_0031))*IT_0038 + (conj(IT_0021)*IT_0027 + IT_0021*conj(IT_0027) + conj
      (IT_0014)*IT_0031 + IT_0014*conj(IT_0031))*IT_0041 + (conj(IT_0014)
      *IT_0027 + IT_0014*conj(IT_0027))*IT_0042) + IT_0043*(conj(IT_0026)
      *IT_0044 + IT_0026*IT_0045 + conj(IT_0035)*IT_0046 + IT_0035*IT_0047) + 
      (conj(IT_0030)*IT_0044 + IT_0030*IT_0045 + conj(IT_0036)*IT_0046 + IT_0036
      *IT_0047)*IT_0050 + IT_0052 + IT_0058 + (conj(IT_0036)*IT_0054 + IT_0036
      *IT_0055 + conj(IT_0030)*IT_0056 + IT_0030*IT_0057)*IT_0059 + (conj
      (IT_0035)*IT_0054 + IT_0035*IT_0055 + conj(IT_0026)*IT_0056 + IT_0026
      *IT_0057)*IT_0062;
    return create_ccomplex_return(IT_0063);
}

