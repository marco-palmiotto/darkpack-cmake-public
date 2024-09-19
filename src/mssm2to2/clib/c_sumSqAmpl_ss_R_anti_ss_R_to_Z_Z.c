#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_anti_ss_R_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_anti_ss_R_to_Z_Z(
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
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_Z, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_ss_R*Gamma_sr + reg_prop, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*IT_0001;
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0005;
    const ccomplex_t IT_0007 = (-0.666666666666667)*IT_0006;
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = IT_0002*IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_ss_R*Gamma_sr + reg_prop, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = IT_0008*IT_0012;
    const ccomplex_t IT_0014 = pow(m_Z, -2);
    const ccomplex_t IT_0015 = s_13*s_23*IT_0014;
    const ccomplex_t IT_0016 = pow(s_24, 2);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = pow(m_ss_R, 2);
    const ccomplex_t IT_0019 = (s_12 + -IT_0015)*(IT_0017 + -IT_0018);
    const ccomplex_t IT_0020 = pow(e_em, 2);
    const ccomplex_t IT_0021 = cpow(IT_0003, -4);
    const ccomplex_t IT_0022 = cpow(IT_0005, 2);
    const ccomplex_t IT_0023 = tan(theta_W);
    const ccomplex_t IT_0024 = cpow(IT_0023, 2);
    const ccomplex_t IT_0025 = cpow(1 + IT_0024, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0020*IT_0021*IT_0022
      *IT_0025;
    const ccomplex_t IT_0027 = 0.222222222222222*IT_0026;
    const ccomplex_t IT_0028 = cos(beta);
    const ccomplex_t IT_0029 = cos(alpha);
    const ccomplex_t IT_0030 = IT_0005*IT_0021*IT_0028*IT_0029;
    const ccomplex_t IT_0031 = sin(beta);
    const ccomplex_t IT_0032 = sin(alpha);
    const ccomplex_t IT_0033 = IT_0005*IT_0021*IT_0031*IT_0032;
    const ccomplex_t IT_0034 = cpow(IT_0005, -3);
    const ccomplex_t IT_0035 = IT_0028*IT_0029*IT_0034;
    const ccomplex_t IT_0036 = IT_0031*IT_0032*IT_0034;
    const ccomplex_t IT_0037 = cpow(IT_0003, -2);
    const ccomplex_t IT_0038 = cpow(IT_0005, -1);
    const ccomplex_t IT_0039 = IT_0028*IT_0029*IT_0037*IT_0038;
    const ccomplex_t IT_0040 = IT_0031*IT_0032*IT_0037*IT_0038;
    const ccomplex_t IT_0041 = cpow(IT_0005, -2);
    const ccomplex_t IT_0042 = cpow(IT_0037 + IT_0041, -1);
    const ccomplex_t IT_0043 = (IT_0030 + IT_0033 + IT_0035 + IT_0036 + 2
      *IT_0039 + 2*IT_0040)*IT_0042;
    const ccomplex_t IT_0044 = m_W*e_em;
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = pow(m_W, -1);
    const ccomplex_t IT_0048 = pow(m_s, 2);
    const ccomplex_t IT_0049 = cpow(IT_0028, -1);
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*e_em*(m_W*IT_0005*(IT_0028
      *IT_0029 + -IT_0031*IT_0032)*IT_0037 + (-3)*IT_0029*IT_0038*IT_0047
      *IT_0048*IT_0049);
    const ccomplex_t IT_0051 = 0.333333333333333*IT_0050;
    const ccomplex_t IT_0052 = IT_0046*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_12 + (-2)*IT_0018 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0005*IT_0021*IT_0029*IT_0031;
    const ccomplex_t IT_0057 = IT_0005*IT_0021*IT_0028*IT_0032;
    const ccomplex_t IT_0058 = IT_0029*IT_0031*IT_0034;
    const ccomplex_t IT_0059 = IT_0028*IT_0032*IT_0034;
    const ccomplex_t IT_0060 = IT_0029*IT_0031*IT_0037*IT_0038;
    const ccomplex_t IT_0061 = IT_0028*IT_0032*IT_0037*IT_0038;
    const ccomplex_t IT_0062 = IT_0042*(IT_0056 + -IT_0057 + IT_0058 + 
      -IT_0059 + 2*IT_0060 + (-2)*IT_0061);
    const ccomplex_t IT_0063 = IT_0044*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*e_em*(m_W*IT_0005*(IT_0029
      *IT_0031 + IT_0028*IT_0032)*IT_0037 + (-3)*IT_0032*IT_0038*IT_0047*IT_0048
      *IT_0049);
    const ccomplex_t IT_0066 = (-0.333333333333333)*IT_0065;
    const ccomplex_t IT_0067 = IT_0064*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0027 + -IT_0055 + -IT_0070;
    const ccomplex_t IT_0072 = pow(m_Z, -4);
    const ccomplex_t IT_0073 = s_13*s_24*s_34*IT_0072;
    const ccomplex_t IT_0074 = s_14*s_24*IT_0014;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = s_12 + -IT_0015 + IT_0073 + IT_0075;
    const ccomplex_t IT_0077 = pow(s_34, 2);
    const ccomplex_t IT_0078 = IT_0072*IT_0077;
    const ccomplex_t IT_0079 = -IT_0013;
    const ccomplex_t IT_0080 = s_34*IT_0014;
    const ccomplex_t IT_0081 = s_23*s_34*IT_0014;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = s_24 + IT_0082;
    const ccomplex_t IT_0084 = IT_0080*IT_0083;
    const ccomplex_t IT_0085 = pow(s_23, 2);
    const ccomplex_t IT_0086 = IT_0014*IT_0085;
    const ccomplex_t IT_0087 = -IT_0018;
    const ccomplex_t IT_0088 = s_23*s_24*s_34*IT_0072;
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = IT_0017 + IT_0086 + IT_0087 + IT_0089;
    const ccomplex_t IT_0091 = pow(s_13, 2);
    const ccomplex_t IT_0092 = IT_0014*IT_0091;
    const ccomplex_t IT_0093 = (IT_0017 + -IT_0018)*(IT_0018 + -IT_0092);
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = s_24*s_34*IT_0014;
    const ccomplex_t IT_0096 = (s_12 + -IT_0015)*(s_23 + -IT_0095);
    const ccomplex_t IT_0097 = 3*IT_0010;
    const ccomplex_t IT_0098 = -(IT_0018 + -IT_0086)*(s_23 + -IT_0095);
    const ccomplex_t IT_0099 = IT_0014*IT_0077;
    const ccomplex_t IT_0100 = -(IT_0018 + -IT_0086)*(IT_0000 + -IT_0099);
    const ccomplex_t IT_0101 = (-3)*conj(IT_0071);
    const ccomplex_t IT_0102 = (IT_0017 + IT_0087)*(IT_0086 + IT_0087);
    const ccomplex_t IT_0103 = (-3)*IT_0010*(conj(IT_0013)*IT_0019 + -conj
      (IT_0071)*IT_0076) + 6*IT_0071*(0.5*conj(IT_0010)*IT_0076 + conj(IT_0071)*
      (1 + 0.5*IT_0078) + (-0.5)*conj(IT_0079)*IT_0084 + (-0.5)*conj(IT_0013)
      *IT_0090) + (conj(IT_0010)*IT_0094 + conj(IT_0079)*IT_0096)*IT_0097 + 3
      *IT_0079*(conj(IT_0010)*IT_0096 + -conj(IT_0013)*IT_0098 + -conj(IT_0079)
      *IT_0100 + 0.333333333333333*IT_0084*IT_0101) + (-3)*IT_0013*(conj(IT_0010
      )*IT_0019 + conj(IT_0079)*IT_0098 + (-0.333333333333333)*IT_0090*IT_0101 +
       -conj(IT_0013)*IT_0102);
    return create_ccomplex_return(IT_0103);
}

