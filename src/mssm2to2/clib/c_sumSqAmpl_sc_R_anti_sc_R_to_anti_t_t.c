#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_anti_sc_R_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_anti_sc_R_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
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
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_sc_R, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -2);
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = pow(m_c, 2);
    const ccomplex_t IT_0011 = cpow(IT_0005, -1);
    const ccomplex_t IT_0012 = cpow(IT_0008, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*(IT_0004
      *IT_0005 + IT_0006*IT_0007)*IT_0008 + (-1.5)*IT_0004*IT_0009*IT_0010
      *IT_0011*IT_0012);
    const ccomplex_t IT_0014 = 0.666666666666667*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_t*e_em*IT_0004*IT_0009
      *IT_0011*IT_0012;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = IT_0001*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*(IT_0004
      *IT_0006 + -IT_0005*IT_0007)*IT_0008 + 1.5*IT_0007*IT_0009*IT_0010*IT_0011
      *IT_0012);
    const ccomplex_t IT_0021 = (-0.666666666666667)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_t*e_em*IT_0007*IT_0009
      *IT_0011*IT_0012;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0019 + -IT_0027;
    const ccomplex_t IT_0029 = pow(m_t, 2);
    const ccomplex_t IT_0030 = cpow(IT_0002, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*e_em*IT_0008*IT_0030;
    const ccomplex_t IT_0032 = 0.666666666666667*IT_0031;
    const ccomplex_t IT_0033 = m_t*IT_0032;
    const ccomplex_t IT_0034 = 1.33333333333333*IT_0031;
    const ccomplex_t IT_0035 = m_t*IT_0034;
    const ccomplex_t IT_0036 = (-0.666666666666667)*IT_0031;
    const ccomplex_t IT_0037 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = IT_0008*IT_0030;
    const ccomplex_t IT_0041 = e_em*IT_0040;
    const ccomplex_t IT_0042 = IT_0002*IT_0012;
    const ccomplex_t IT_0043 = e_em*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0041 + (-3)*IT_0043);
    const ccomplex_t IT_0045 = (-0.166666666666667)*IT_0044;
    const ccomplex_t IT_0046 = IT_0037*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = (IT_0033 + -IT_0035)*(IT_0039 + -IT_0047);
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = -(IT_0033 + -IT_0035)*(IT_0039 + -IT_0047);
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = (-18)*IT_0029;
    const ccomplex_t IT_0055 = 18*s_34;
    const ccomplex_t IT_0056 = IT_0054 + IT_0055;
    const ccomplex_t IT_0057 = IT_0034*IT_0047;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = m_t*s_13;
    const ccomplex_t IT_0060 = s_34*IT_0000;
    const ccomplex_t IT_0061 = (-18)*IT_0060;
    const ccomplex_t IT_0062 = s_13*s_14;
    const ccomplex_t IT_0063 = 36*IT_0062;
    const ccomplex_t IT_0064 = IT_0061 + IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0066 = 0.666666666666667*IT_0065;
    const ccomplex_t IT_0067 = (-1.33333333333333)*IT_0065;
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = cpow(s_12 + IT_0000 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0073 = (-2)*IT_0072;
    const ccomplex_t IT_0074 = IT_0069*IT_0072;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0073*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0071 + (-0.0833333333333333)*IT_0076;
    const ccomplex_t IT_0078 = cpow(IT_0077, 2);
    const ccomplex_t IT_0079 = (-36)*IT_0060;
    const ccomplex_t IT_0080 = IT_0000*IT_0029;
    const ccomplex_t IT_0081 = (-36)*IT_0080;
    const ccomplex_t IT_0082 = 72*IT_0062;
    const ccomplex_t IT_0083 = IT_0079 + IT_0081 + IT_0082;
    const ccomplex_t IT_0084 = 0.25*IT_0076;
    const ccomplex_t IT_0085 = IT_0034*IT_0039;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = (-6)*IT_0060;
    const ccomplex_t IT_0088 = (-6)*IT_0080;
    const ccomplex_t IT_0089 = 12*IT_0062;
    const ccomplex_t IT_0090 = IT_0087 + IT_0088 + IT_0089;
    const ccomplex_t IT_0091 = -s_14;
    const ccomplex_t IT_0092 = s_13 + IT_0091;
    const ccomplex_t IT_0093 = m_t*IT_0092;
    const ccomplex_t IT_0094 = 18*IT_0028;
    const ccomplex_t IT_0095 = IT_0058 + IT_0086;
    const ccomplex_t IT_0096 = 18*conj(IT_0028);
    const ccomplex_t IT_0097 = (-18)*IT_0050;
    const ccomplex_t IT_0098 = 2*(IT_0028 + -conj(IT_0028) + 1./2*IT_0050 + 
      -1./2*conj(IT_0050) + 1./2*IT_0053 + -1./2*conj(IT_0053))*IT_0084 + (
      -0.166666666666667)*(conj(IT_0058) + conj(IT_0086))*IT_0094 + (
      -0.166666666666667)*IT_0095*IT_0096 + (-0.166666666666667)*IT_0077*((-36)
      *IT_0028 + 36*conj(IT_0028) + 18*conj(IT_0050) + (-18)*IT_0053 + 18*conj
      (IT_0053) + IT_0097);
    const ccomplex_t IT_0099 = (-18)*IT_0080;
    const ccomplex_t IT_0100 = IT_0061 + IT_0063 + IT_0099;
    const ccomplex_t IT_0101 = m_t*s_14;
    const ccomplex_t IT_0102 = (-18)*IT_0058;
    const ccomplex_t IT_0103 = (-18)*conj(IT_0050);
    const ccomplex_t IT_0104 = (-18)*conj(IT_0058);
    const ccomplex_t IT_0105 = -IT_0077;
    const ccomplex_t IT_0106 = IT_0028*IT_0056;
    const ccomplex_t IT_0107 = conj(IT_0028)*(IT_0028*(36*s_34 + (-36)*IT_0029
      ) + (IT_0050 + IT_0053)*IT_0056) + conj(IT_0058)*(18*IT_0050*IT_0059 +
       IT_0058*IT_0064) + -IT_0078*IT_0083 + (-2)*IT_0084*(IT_0077*((-12)
      *IT_0060 + 24*IT_0062 + (-12)*IT_0080) + 0.5*IT_0083*IT_0084 + 0.5*
      (IT_0058 + -conj(IT_0058) + IT_0086 + -conj(IT_0086))*IT_0090) + (-6)
      *IT_0093*IT_0098 + conj(IT_0086)*(18*IT_0053*IT_0059 + IT_0064*IT_0086 +
       IT_0077*IT_0100 + IT_0097*IT_0101 + IT_0080*IT_0102) + IT_0086*(IT_0101
      *IT_0103 + IT_0080*IT_0104) + IT_0053*(IT_0029*IT_0103 + IT_0101*IT_0104) 
      + IT_0100*(conj(IT_0058)*IT_0077 + IT_0095*IT_0105) + conj(IT_0050)*
      (IT_0050*IT_0055 + 18*IT_0058*IT_0059 + IT_0106) + conj(IT_0053)*(IT_0053
      *IT_0055 + 18*IT_0059*IT_0086 + IT_0029*IT_0097 + IT_0101*IT_0102 +
       IT_0106);
    return create_ccomplex_return(IT_0107);
}

