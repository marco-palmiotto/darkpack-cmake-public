#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_L_to_anti_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_L_to_anti_mu_mu(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*m_mu*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = pow(m_s, 2);
    const ccomplex_t IT_0009 = cos(alpha);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0003*IT_0005*IT_0008 + (-0.5)*m_W*(IT_0001*IT_0003 + IT_0009*IT_0010)*
      (IT_0005 + 1./3*IT_0004*IT_0012));
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = pow(m_ss_L, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*m_mu*IT_0000*IT_0002
      *IT_0005*IT_0009;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0005*IT_0008*IT_0009 + (-0.5)*m_W*(IT_0001*IT_0009 + -IT_0003*IT_0010)
      *(IT_0005 + 1./3*IT_0004*IT_0012));
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + (-2)*IT_0015 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0018 + -IT_0026;
    const ccomplex_t IT_0028 = pow(m_mu, 2);
    const ccomplex_t IT_0029 = cpow(IT_0011, -1);
    const ccomplex_t IT_0030 = (0 + _Complex_I*3)*e_em*(IT_0005*IT_0011 + 1./3
      *IT_0004*IT_0029);
    const ccomplex_t IT_0031 = 0.333333333333333*IT_0030;
    const ccomplex_t IT_0032 = IT_0004*IT_0029;
    const ccomplex_t IT_0033 = e_em*IT_0032;
    const ccomplex_t IT_0034 = IT_0005*IT_0011;
    const ccomplex_t IT_0035 = e_em*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0033 + -IT_0035);
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0031*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0029;
    const ccomplex_t IT_0044 = IT_0038*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0031*IT_0045;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = -s_14;
    const ccomplex_t IT_0049 = s_13 + IT_0048;
    const ccomplex_t IT_0050 = m_mu*IT_0049;
    const ccomplex_t IT_0051 = 0.166666666666667*IT_0030;
    const ccomplex_t IT_0052 = m_mu*IT_0051;
    const ccomplex_t IT_0053 = m_mu*IT_0031;
    const ccomplex_t IT_0054 = -(IT_0040 + -IT_0045)*(IT_0052 + -IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = -IT_0028;
    const ccomplex_t IT_0058 = s_34 + IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = 0.666666666666667*IT_0059;
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = cpow(s_12 + IT_0015 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = s_13*s_14;
    const ccomplex_t IT_0068 = IT_0015*IT_0028;
    const ccomplex_t IT_0069 = s_34*IT_0015;
    const ccomplex_t IT_0070 = (-6)*IT_0069;
    const ccomplex_t IT_0071 = (-6)*IT_0068;
    const ccomplex_t IT_0072 = 12*IT_0067;
    const ccomplex_t IT_0073 = IT_0070 + IT_0071 + IT_0072;
    const ccomplex_t IT_0074 = (-6)*IT_0056;
    const ccomplex_t IT_0075 = IT_0066*(24*IT_0067 + (-12)*IT_0068 + (-12)
      *IT_0069) + (IT_0042 + -conj(IT_0042) + IT_0047 + -conj(IT_0047))*IT_0073 
      + -IT_0050*((-12)*IT_0027 + 12*conj(IT_0027) + (-6)*IT_0055 + 6*conj
      (IT_0055) + 6*conj(IT_0056) + IT_0074);
    const ccomplex_t IT_0076 = m_mu*s_14;
    const ccomplex_t IT_0077 = (-6)*IT_0042;
    const ccomplex_t IT_0078 = m_mu*s_13;
    const ccomplex_t IT_0079 = IT_0070 + IT_0072;
    const ccomplex_t IT_0080 = conj(IT_0027)*IT_0050;
    const ccomplex_t IT_0081 = (-6)*conj(IT_0042);
    const ccomplex_t IT_0082 = (-6)*conj(IT_0056);
    const ccomplex_t IT_0083 = 0.166666666666667*IT_0076;
    const ccomplex_t IT_0084 = 6*s_34;
    const ccomplex_t IT_0085 = conj(IT_0027)*IT_0058;
    const ccomplex_t IT_0086 = IT_0027*(conj(IT_0027)*(12*s_34 + (-12)*IT_0028
      ) + 6*(conj(IT_0042) + conj(IT_0047))*IT_0050 + 6*(conj(IT_0055) + conj
      (IT_0056))*IT_0058) + -IT_0066*IT_0075 + conj(IT_0047)*(IT_0074*IT_0076 +
       IT_0068*IT_0077) + conj(IT_0055)*(IT_0028*IT_0074 + IT_0076*IT_0077) + 6
      *IT_0042*(conj(IT_0056)*IT_0078 + 0.166666666666667*conj(IT_0042)*IT_0079 
      + IT_0080) + 6*IT_0047*(conj(IT_0055)*IT_0078 + 0.166666666666667*conj
      (IT_0047)*IT_0079 + IT_0080 + 0.166666666666667*IT_0068*IT_0081 + IT_0082
      *IT_0083) + 6*IT_0055*(conj(IT_0047)*IT_0078 + 0.166666666666667*IT_0028
      *IT_0082 + IT_0081*IT_0083 + 0.166666666666667*conj(IT_0055)*IT_0084 +
       IT_0085) + 6*IT_0056*(conj(IT_0042)*IT_0078 + 0.166666666666667*conj
      (IT_0056)*IT_0084 + IT_0085);
    return create_ccomplex_return(IT_0086);
}

