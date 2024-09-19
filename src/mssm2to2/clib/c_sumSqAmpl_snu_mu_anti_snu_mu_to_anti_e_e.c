#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_mu_anti_snu_mu_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_mu_anti_snu_mu_to_anti_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0003;
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0002;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0003;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0007 + -IT_0009);
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0005*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = (-0.5)*IT_0004;
    const ccomplex_t IT_0019 = m_e*IT_0018;
    const ccomplex_t IT_0020 = m_e*IT_0005;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0022 = IT_0013*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = (IT_0019 + -IT_0020)*(IT_0015 + -IT_0023);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = m_e*s_13;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = cos(beta);
    const ccomplex_t IT_0029 = cpow(IT_0028, -1);
    const ccomplex_t IT_0030 = sin(alpha);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*m_e*e_em*IT_0003*IT_0027
      *IT_0029*IT_0030;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = cpow(IT_0000, -2);
    const ccomplex_t IT_0034 = IT_0002*IT_0033;
    const ccomplex_t IT_0035 = IT_0003 + IT_0034;
    const ccomplex_t IT_0036 = cos(alpha);
    const ccomplex_t IT_0037 = sin(beta);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_W*e_em*IT_0035*(IT_0028
      *IT_0030 + IT_0036*IT_0037);
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0032*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*m_e*e_em*IT_0003*IT_0027
      *IT_0029*IT_0036;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*m_W*e_em*IT_0035*(IT_0028
      *IT_0036 + -IT_0030*IT_0037);
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + (-2)*IT_0012 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = -IT_0043 + -IT_0051;
    const ccomplex_t IT_0053 = -s_14;
    const ccomplex_t IT_0054 = s_13 + IT_0053;
    const ccomplex_t IT_0055 = m_e*IT_0054;
    const ccomplex_t IT_0056 = s_13*s_14;
    const ccomplex_t IT_0057 = s_34*IT_0012;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = pow(m_e, 2);
    const ccomplex_t IT_0061 = IT_0005*IT_0023;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = -IT_0060;
    const ccomplex_t IT_0064 = s_34 + IT_0063;
    const ccomplex_t IT_0065 = 2*s_34;
    const ccomplex_t IT_0066 = m_e*s_14;
    const ccomplex_t IT_0067 = (-2)*IT_0017;
    const ccomplex_t IT_0068 = IT_0012*IT_0060;
    const ccomplex_t IT_0069 = (-2)*conj(IT_0025);
    const ccomplex_t IT_0070 = (-2)*conj(IT_0017);
    const ccomplex_t IT_0071 = 0.5*IT_0025;
    const ccomplex_t IT_0072 = 2*conj(IT_0017)*(IT_0025*IT_0026 + IT_0052
      *IT_0055 + 2*IT_0017*IT_0059) + (-2)*conj(IT_0024)*(IT_0025*IT_0060 + 
      -IT_0026*IT_0062 + -IT_0052*IT_0064 + (-0.5)*IT_0024*IT_0065 + (-0.5)
      *IT_0066*IT_0067) + 2*conj(IT_0062)*(IT_0024*IT_0026 + IT_0052*IT_0055 + 2
      *IT_0059*IT_0062 + -IT_0025*IT_0066 + 0.5*IT_0067*IT_0068) + IT_0024*
      (IT_0060*IT_0069 + IT_0066*IT_0070) + IT_0062*(IT_0066*IT_0069 + IT_0068
      *IT_0070) + 2*conj(IT_0052)*(IT_0055*(IT_0017 + IT_0062) + 2*IT_0064*(0.5
      *IT_0024 + IT_0052 + IT_0071)) + 2*conj(IT_0025)*(IT_0017*IT_0026 +
       IT_0052*IT_0064 + IT_0065*IT_0071);
    return create_ccomplex_return(IT_0072);
}

