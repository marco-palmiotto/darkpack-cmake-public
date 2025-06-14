#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_anti_sd_R_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_anti_sd_R_to_anti_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
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
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = pow(m_d, 2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0001*IT_0003 + -IT_0010*IT_0011) + (-3)*IT_0000*IT_0002*IT_0003
      *IT_0005*IT_0012);
    const ccomplex_t IT_0014 = 0.333333333333333*IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_sd_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0005*IT_0011;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0003*IT_0010 + IT_0001*IT_0011) + (-3)*IT_0000*IT_0002*IT_0005*IT_0011
      *IT_0012);
    const ccomplex_t IT_0023 = (-0.333333333333333)*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0019 + -IT_0027;
    const ccomplex_t IT_0029 = pow(m_e, 2);
    const ccomplex_t IT_0030 = cpow(IT_0008, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0030;
    const ccomplex_t IT_0032 = (-0.666666666666667)*IT_0031;
    const ccomplex_t IT_0033 = IT_0004*IT_0030;
    const ccomplex_t IT_0034 = e_em*IT_0033;
    const ccomplex_t IT_0035 = IT_0005*IT_0008;
    const ccomplex_t IT_0036 = e_em*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0034 + -IT_0036);
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0032*IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0031*IT_0039;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0032*IT_0045;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = -s_14;
    const ccomplex_t IT_0049 = s_13 + IT_0048;
    const ccomplex_t IT_0050 = m_e*IT_0049;
    const ccomplex_t IT_0051 = m_e*IT_0032;
    const ccomplex_t IT_0052 = (-0.333333333333333)*IT_0031;
    const ccomplex_t IT_0053 = m_e*IT_0052;
    const ccomplex_t IT_0054 = (IT_0041 + -IT_0045)*(IT_0051 + -IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = -IT_0029;
    const ccomplex_t IT_0057 = s_34 + IT_0056;
    const ccomplex_t IT_0058 = m_e*s_14;
    const ccomplex_t IT_0059 = (-6)*IT_0055;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = 0.666666666666667*IT_0060;
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = s_13*s_14;
    const ccomplex_t IT_0069 = IT_0016*IT_0029;
    const ccomplex_t IT_0070 = s_34*IT_0016;
    const ccomplex_t IT_0071 = (-6)*IT_0070;
    const ccomplex_t IT_0072 = (-6)*IT_0069;
    const ccomplex_t IT_0073 = 12*IT_0068;
    const ccomplex_t IT_0074 = IT_0071 + IT_0072 + IT_0073;
    const ccomplex_t IT_0075 = -conj(IT_0047);
    const ccomplex_t IT_0076 = -IT_0050*((-12)*IT_0028 + 12*conj(IT_0028) + (
      -6)*IT_0054 + 6*conj(IT_0054) + 6*conj(IT_0055) + IT_0059) + IT_0067*(24
      *IT_0068 + (-12)*IT_0069 + (-12)*IT_0070) + IT_0074*(IT_0043 + -conj
      (IT_0043) + IT_0047 + IT_0075);
    const ccomplex_t IT_0077 = m_e*s_13;
    const ccomplex_t IT_0078 = IT_0071 + IT_0073;
    const ccomplex_t IT_0079 = conj(IT_0028)*IT_0050;
    const ccomplex_t IT_0080 = (-6)*conj(IT_0055);
    const ccomplex_t IT_0081 = 6*s_34;
    const ccomplex_t IT_0082 = conj(IT_0028)*IT_0057;
    const ccomplex_t IT_0083 = IT_0028*(conj(IT_0028)*(12*s_34 + (-12)*IT_0029
      ) + 6*(conj(IT_0043) + conj(IT_0047))*IT_0050 + 6*(conj(IT_0054) + conj
      (IT_0055))*IT_0057) + (IT_0029*conj(IT_0054) + conj(IT_0043)*IT_0058)
      *IT_0059 + -IT_0067*IT_0076 + (-6)*IT_0047*(conj(IT_0054)*IT_0058 + conj
      (IT_0043)*IT_0069 + -conj(IT_0055)*IT_0077 + (-0.166666666666667)*conj
      (IT_0047)*IT_0078 + -IT_0079) + 6*IT_0043*(IT_0069*IT_0075 + conj(IT_0054)
      *IT_0077 + 0.166666666666667*conj(IT_0043)*IT_0078 + IT_0079 +
       0.166666666666667*IT_0058*IT_0080) + 6*IT_0054*(IT_0058*IT_0075 + conj
      (IT_0043)*IT_0077 + 0.166666666666667*IT_0029*IT_0080 + 0.166666666666667
      *conj(IT_0054)*IT_0081 + IT_0082) + 6*IT_0055*(conj(IT_0047)*IT_0077 +
       0.166666666666667*conj(IT_0055)*IT_0081 + IT_0082);
    return create_ccomplex_return(IT_0083);
}

