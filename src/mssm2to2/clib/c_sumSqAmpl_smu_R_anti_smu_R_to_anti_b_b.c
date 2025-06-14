#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_smu_R_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_smu_R_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = pow(m_mu, 2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0001*IT_0003 + -IT_0010*IT_0011) + -IT_0000*IT_0002*IT_0003*IT_0005
      *IT_0012);
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = pow(m_smu_R, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + (-2)*IT_0015 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0005*IT_0011;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0003*IT_0010 + IT_0001*IT_0011) + -IT_0000*IT_0002*IT_0005*IT_0011
      *IT_0012);
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0018 + -IT_0026;
    const ccomplex_t IT_0028 = pow(m_b, 2);
    const ccomplex_t IT_0029 = cpow(IT_0008, -1);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0029;
    const ccomplex_t IT_0031 = (-2)*IT_0030;
    const ccomplex_t IT_0032 = IT_0004*IT_0029;
    const ccomplex_t IT_0033 = e_em*IT_0032;
    const ccomplex_t IT_0034 = IT_0005*IT_0008;
    const ccomplex_t IT_0035 = e_em*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0033 + 3*IT_0035);
    const ccomplex_t IT_0037 = (-0.166666666666667)*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0031*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = 0.333333333333333*IT_0030;
    const ccomplex_t IT_0044 = IT_0038*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0031*IT_0045;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = -s_14;
    const ccomplex_t IT_0049 = s_13 + IT_0048;
    const ccomplex_t IT_0050 = m_b*IT_0049;
    const ccomplex_t IT_0051 = -IT_0030;
    const ccomplex_t IT_0052 = m_b*IT_0051;
    const ccomplex_t IT_0053 = m_b*IT_0031;
    const ccomplex_t IT_0054 = -(IT_0040 + -IT_0045)*(IT_0052 + -IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = (IT_0040 + -IT_0045)*(IT_0052 + -IT_0053);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = -IT_0028;
    const ccomplex_t IT_0061 = s_34 + IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0063 = 2*IT_0062;
    const ccomplex_t IT_0064 = (-0.333333333333333)*IT_0062;
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = cpow(s_12 + IT_0015 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = s_13*s_14;
    const ccomplex_t IT_0071 = IT_0015*IT_0028;
    const ccomplex_t IT_0072 = s_34*IT_0015;
    const ccomplex_t IT_0073 = (-6)*IT_0072;
    const ccomplex_t IT_0074 = (-6)*IT_0071;
    const ccomplex_t IT_0075 = 12*IT_0070;
    const ccomplex_t IT_0076 = IT_0073 + IT_0074 + IT_0075;
    const ccomplex_t IT_0077 = (-6)*IT_0056;
    const ccomplex_t IT_0078 = IT_0069*(24*IT_0070 + (-12)*IT_0071 + (-12)
      *IT_0072) + (IT_0042 + -conj(IT_0042) + IT_0047 + -conj(IT_0047))*IT_0076 
      + -IT_0050*((-12)*IT_0027 + 12*conj(IT_0027) + 6*conj(IT_0056) + (-6)
      *IT_0059 + 6*conj(IT_0059) + IT_0077);
    const ccomplex_t IT_0079 = m_b*s_14;
    const ccomplex_t IT_0080 = (-6)*IT_0042;
    const ccomplex_t IT_0081 = m_b*s_13;
    const ccomplex_t IT_0082 = IT_0073 + IT_0075;
    const ccomplex_t IT_0083 = conj(IT_0027)*IT_0050;
    const ccomplex_t IT_0084 = (-6)*conj(IT_0042);
    const ccomplex_t IT_0085 = (-6)*conj(IT_0056);
    const ccomplex_t IT_0086 = 0.166666666666667*IT_0079;
    const ccomplex_t IT_0087 = 6*s_34;
    const ccomplex_t IT_0088 = conj(IT_0027)*IT_0061;
    const ccomplex_t IT_0089 = IT_0027*(conj(IT_0027)*(12*s_34 + (-12)*IT_0028
      ) + 6*(conj(IT_0042) + conj(IT_0047))*IT_0050 + 6*(conj(IT_0056) + conj
      (IT_0059))*IT_0061) + -IT_0069*IT_0078 + conj(IT_0047)*(IT_0077*IT_0079 +
       IT_0071*IT_0080) + conj(IT_0059)*(IT_0028*IT_0077 + IT_0079*IT_0080) + 6
      *IT_0042*(conj(IT_0056)*IT_0081 + 0.166666666666667*conj(IT_0042)*IT_0082 
      + IT_0083) + 6*IT_0047*(conj(IT_0059)*IT_0081 + 0.166666666666667*conj
      (IT_0047)*IT_0082 + IT_0083 + 0.166666666666667*IT_0071*IT_0084 + IT_0085
      *IT_0086) + 6*IT_0056*(conj(IT_0042)*IT_0081 + 0.166666666666667*conj
      (IT_0056)*IT_0087 + IT_0088) + 6*IT_0059*(conj(IT_0047)*IT_0081 +
       0.166666666666667*IT_0028*IT_0085 + IT_0084*IT_0086 + 0.166666666666667
      *conj(IT_0059)*IT_0087 + IT_0088);
    return create_ccomplex_return(IT_0089);
}

