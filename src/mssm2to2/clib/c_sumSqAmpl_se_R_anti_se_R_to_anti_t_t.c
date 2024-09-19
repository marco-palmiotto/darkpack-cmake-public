#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_se_R_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_se_R_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -2);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = pow(m_e, 2);
    const ccomplex_t IT_0009 = cpow(IT_0004, -1);
    const ccomplex_t IT_0010 = cpow(IT_0006, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*(m_W*IT_0001*(IT_0002
      *IT_0003 + IT_0004*IT_0005)*IT_0006 + -IT_0005*IT_0007*IT_0008*IT_0009
      *IT_0010);
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0003, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_t*e_em*IT_0002*IT_0007
      *IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = pow(m_se_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(m_W*IT_0001*(IT_0002
      *IT_0004 + -IT_0003*IT_0005)*IT_0006 + -IT_0002*IT_0007*IT_0008*IT_0009
      *IT_0010);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_t*e_em*IT_0005*IT_0007
      *IT_0010*IT_0013;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0020 + -IT_0027;
    const ccomplex_t IT_0029 = pow(m_t, 2);
    const ccomplex_t IT_0030 = cpow(IT_0000, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0030;
    const ccomplex_t IT_0032 = (-2)*IT_0031;
    const ccomplex_t IT_0033 = IT_0006*IT_0030;
    const ccomplex_t IT_0034 = e_em*IT_0033;
    const ccomplex_t IT_0035 = IT_0000*IT_0010;
    const ccomplex_t IT_0036 = e_em*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0034 + (-3)*IT_0036);
    const ccomplex_t IT_0038 = (-0.166666666666667)*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0032*IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = (-0.666666666666667)*IT_0031;
    const ccomplex_t IT_0045 = IT_0039*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0032*IT_0046;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = -s_14;
    const ccomplex_t IT_0050 = s_13 + IT_0049;
    const ccomplex_t IT_0051 = m_t*IT_0050;
    const ccomplex_t IT_0052 = -IT_0031;
    const ccomplex_t IT_0053 = m_t*IT_0052;
    const ccomplex_t IT_0054 = m_t*IT_0032;
    const ccomplex_t IT_0055 = -(IT_0041 + -IT_0046)*(IT_0053 + -IT_0054);
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = (IT_0041 + -IT_0046)*(IT_0053 + -IT_0054);
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = -IT_0029;
    const ccomplex_t IT_0062 = s_34 + IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0064 = 2*IT_0063;
    const ccomplex_t IT_0065 = 0.666666666666667*IT_0063;
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = cpow(s_12 + IT_0017 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = s_13*s_14;
    const ccomplex_t IT_0072 = IT_0017*IT_0029;
    const ccomplex_t IT_0073 = s_34*IT_0017;
    const ccomplex_t IT_0074 = (-6)*IT_0073;
    const ccomplex_t IT_0075 = (-6)*IT_0072;
    const ccomplex_t IT_0076 = 12*IT_0071;
    const ccomplex_t IT_0077 = IT_0074 + IT_0075 + IT_0076;
    const ccomplex_t IT_0078 = (-6)*IT_0057;
    const ccomplex_t IT_0079 = IT_0070*(24*IT_0071 + (-12)*IT_0072 + (-12)
      *IT_0073) + (IT_0043 + -conj(IT_0043) + IT_0048 + -conj(IT_0048))*IT_0077 
      + -IT_0051*((-12)*IT_0028 + 12*conj(IT_0028) + 6*conj(IT_0057) + (-6)
      *IT_0060 + 6*conj(IT_0060) + IT_0078);
    const ccomplex_t IT_0080 = m_t*s_14;
    const ccomplex_t IT_0081 = (-6)*IT_0043;
    const ccomplex_t IT_0082 = m_t*s_13;
    const ccomplex_t IT_0083 = IT_0074 + IT_0076;
    const ccomplex_t IT_0084 = conj(IT_0028)*IT_0051;
    const ccomplex_t IT_0085 = (-6)*conj(IT_0043);
    const ccomplex_t IT_0086 = (-6)*conj(IT_0057);
    const ccomplex_t IT_0087 = 0.166666666666667*IT_0080;
    const ccomplex_t IT_0088 = 6*s_34;
    const ccomplex_t IT_0089 = conj(IT_0028)*IT_0062;
    const ccomplex_t IT_0090 = IT_0028*(conj(IT_0028)*(12*s_34 + (-12)*IT_0029
      ) + 6*(conj(IT_0043) + conj(IT_0048))*IT_0051 + 6*(conj(IT_0057) + conj
      (IT_0060))*IT_0062) + -IT_0070*IT_0079 + conj(IT_0048)*(IT_0078*IT_0080 +
       IT_0072*IT_0081) + conj(IT_0060)*(IT_0029*IT_0078 + IT_0080*IT_0081) + 6
      *IT_0043*(conj(IT_0057)*IT_0082 + 0.166666666666667*conj(IT_0043)*IT_0083 
      + IT_0084) + 6*IT_0048*(conj(IT_0060)*IT_0082 + 0.166666666666667*conj
      (IT_0048)*IT_0083 + IT_0084 + 0.166666666666667*IT_0072*IT_0085 + IT_0086
      *IT_0087) + 6*IT_0057*(conj(IT_0043)*IT_0082 + 0.166666666666667*conj
      (IT_0057)*IT_0088 + IT_0089) + 6*IT_0060*(conj(IT_0048)*IT_0082 +
       0.166666666666667*IT_0029*IT_0086 + IT_0085*IT_0087 + 0.166666666666667
      *conj(IT_0060)*IT_0088 + IT_0089);
    return create_ccomplex_return(IT_0090);
}

