#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_L_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_L_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 3
      *IT_0000*IT_0003);
    const ccomplex_t IT_0005 = 0.166666666666667*IT_0004;
    const ccomplex_t IT_0006 = m_b*IT_0005;
    const ccomplex_t IT_0007 = 0.333333333333333*IT_0004;
    const ccomplex_t IT_0008 = m_b*IT_0007;
    const ccomplex_t IT_0009 = IT_0001*IT_0002;
    const ccomplex_t IT_0010 = e_em*IT_0009;
    const ccomplex_t IT_0011 = IT_0000*IT_0003;
    const ccomplex_t IT_0012 = e_em*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0010 + 3*IT_0012);
    const ccomplex_t IT_0014 = (-0.166666666666667)*IT_0013;
    const ccomplex_t IT_0015 = pow(m_ss_L, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0020 = 0.333333333333333*IT_0019;
    const ccomplex_t IT_0021 = IT_0016*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (IT_0006 + -IT_0008)*(IT_0018 + -IT_0022);
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = pow(m_W, -1);
    const ccomplex_t IT_0026 = cos(beta);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = cos(alpha);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_b*e_em*IT_0003*IT_0025
      *IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = pow(m_s, 2);
    const ccomplex_t IT_0032 = sin(beta);
    const ccomplex_t IT_0033 = sin(alpha);
    const ccomplex_t IT_0034 = cpow(IT_0000, -2);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0025
      *IT_0027*IT_0028*IT_0031 + (-0.5)*m_W*(IT_0026*IT_0028 + -IT_0032*IT_0033)
      *(IT_0003 + 1./3*IT_0002*IT_0034));
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = IT_0030*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_12 + (-2)*IT_0015 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_b*e_em*IT_0003*IT_0025
      *IT_0027*IT_0033;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0025
      *IT_0027*IT_0031*IT_0033 + (-0.5)*m_W*(IT_0028*IT_0032 + IT_0026*IT_0033)*
      (IT_0003 + 1./3*IT_0002*IT_0034));
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = -IT_0040 + -IT_0047;
    const ccomplex_t IT_0049 = pow(m_b, 2);
    const ccomplex_t IT_0050 = (-18)*IT_0049;
    const ccomplex_t IT_0051 = 18*s_34;
    const ccomplex_t IT_0052 = IT_0050 + IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0054 = 0.666666666666667*IT_0053;
    const ccomplex_t IT_0055 = (-0.333333333333333)*IT_0053;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = cpow(s_12 + IT_0015 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0061 = (-2)*IT_0060;
    const ccomplex_t IT_0062 = IT_0057*IT_0060;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0061*IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0059 + (-0.0833333333333333)*IT_0064;
    const ccomplex_t IT_0066 = cpow(IT_0065, 2);
    const ccomplex_t IT_0067 = s_34*IT_0015;
    const ccomplex_t IT_0068 = (-36)*IT_0067;
    const ccomplex_t IT_0069 = IT_0015*IT_0049;
    const ccomplex_t IT_0070 = (-36)*IT_0069;
    const ccomplex_t IT_0071 = s_13*s_14;
    const ccomplex_t IT_0072 = 72*IT_0071;
    const ccomplex_t IT_0073 = IT_0068 + IT_0070 + IT_0072;
    const ccomplex_t IT_0074 = 0.25*IT_0064;
    const ccomplex_t IT_0075 = IT_0007*IT_0018;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0007*IT_0022;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = (-6)*IT_0067;
    const ccomplex_t IT_0080 = (-6)*IT_0069;
    const ccomplex_t IT_0081 = 12*IT_0071;
    const ccomplex_t IT_0082 = IT_0079 + IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = m_b*s_13;
    const ccomplex_t IT_0084 = (-18)*IT_0067;
    const ccomplex_t IT_0085 = 36*IT_0071;
    const ccomplex_t IT_0086 = IT_0084 + IT_0085;
    const ccomplex_t IT_0087 = 18*conj(IT_0024);
    const ccomplex_t IT_0088 = -s_14;
    const ccomplex_t IT_0089 = s_13 + IT_0088;
    const ccomplex_t IT_0090 = m_b*IT_0089;
    const ccomplex_t IT_0091 = 18*IT_0048;
    const ccomplex_t IT_0092 = IT_0076 + IT_0078;
    const ccomplex_t IT_0093 = 18*conj(IT_0048);
    const ccomplex_t IT_0094 = (-18)*IT_0024;
    const ccomplex_t IT_0095 = (IT_0023 + -conj(IT_0023) + IT_0024 + -conj
      (IT_0024) + 2*IT_0048 + (-2)*conj(IT_0048))*IT_0074 + (-0.166666666666667)
      *(conj(IT_0076) + conj(IT_0078))*IT_0091 + (-0.166666666666667)*IT_0092
      *IT_0093 + (-0.166666666666667)*IT_0065*((-18)*IT_0023 + 18*conj(IT_0023) 
      + (-36)*IT_0048 + 36*conj(IT_0048) + IT_0087 + IT_0094);
    const ccomplex_t IT_0096 = m_b*s_14;
    const ccomplex_t IT_0097 = (-18)*IT_0076;
    const ccomplex_t IT_0098 = (-18)*IT_0069;
    const ccomplex_t IT_0099 = IT_0084 + IT_0085 + IT_0098;
    const ccomplex_t IT_0100 = (-18)*conj(IT_0024);
    const ccomplex_t IT_0101 = (-18)*conj(IT_0076);
    const ccomplex_t IT_0102 = -IT_0065;
    const ccomplex_t IT_0103 = conj(IT_0024)*IT_0048*IT_0052 + conj(IT_0048)*
      (IT_0048*(36*s_34 + (-36)*IT_0049) + (IT_0023 + IT_0024)*IT_0052) + 
      -IT_0066*IT_0073 + (-2)*IT_0074*(IT_0065*((-12)*IT_0067 + (-12)*IT_0069 +
       24*IT_0071) + 0.5*IT_0073*IT_0074 + 0.5*(IT_0076 + -conj(IT_0076) +
       IT_0078 + -conj(IT_0078))*IT_0082) + IT_0024*(conj(IT_0024)*IT_0051 + 18
      *conj(IT_0076)*IT_0083) + IT_0076*(conj(IT_0076)*IT_0086 + IT_0083*IT_0087
      ) + (-6)*IT_0090*IT_0095 + conj(IT_0023)*(IT_0023*IT_0051 + IT_0048
      *IT_0052 + 18*IT_0078*IT_0083 + IT_0049*IT_0094 + IT_0096*IT_0097) + conj
      (IT_0078)*(18*IT_0023*IT_0083 + IT_0078*IT_0086 + IT_0094*IT_0096 +
       IT_0069*IT_0097 + IT_0065*IT_0099) + IT_0078*(IT_0096*IT_0100 + IT_0069
      *IT_0101) + IT_0023*(IT_0049*IT_0100 + IT_0096*IT_0101) + IT_0099*(IT_0065
      *conj(IT_0076) + IT_0092*IT_0102);
    return create_ccomplex_return(IT_0103);
}

