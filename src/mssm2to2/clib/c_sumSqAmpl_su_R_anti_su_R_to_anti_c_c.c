#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_su_R_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_su_R_to_anti_c_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = pow(m_u, 2);
    const ccomplex_t IT_0009 = cpow(IT_0004, -1);
    const ccomplex_t IT_0010 = cpow(IT_0006, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*(m_W*IT_0001*(IT_0002
      *IT_0003 + -IT_0004*IT_0005)*IT_0006 + 1.5*IT_0005*IT_0007*IT_0008*IT_0009
      *IT_0010);
    const ccomplex_t IT_0012 = (-0.666666666666667)*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*m_c*e_em*IT_0005*IT_0007
      *IT_0009*IT_0010;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = pow(m_su_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*e_em*(m_W*IT_0001*(IT_0003
      *IT_0004 + IT_0002*IT_0005)*IT_0006 + (-1.5)*IT_0003*IT_0007*IT_0008
      *IT_0009*IT_0010);
    const ccomplex_t IT_0021 = 0.666666666666667*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_c*e_em*IT_0003*IT_0007
      *IT_0009*IT_0010;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0019 + -IT_0027;
    const ccomplex_t IT_0029 = cpow(IT_0000, -1);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0029;
    const ccomplex_t IT_0031 = 0.666666666666667*IT_0030;
    const ccomplex_t IT_0032 = m_c*IT_0031;
    const ccomplex_t IT_0033 = 1.33333333333333*IT_0030;
    const ccomplex_t IT_0034 = m_c*IT_0033;
    const ccomplex_t IT_0035 = (-0.666666666666667)*IT_0030;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0006*IT_0029;
    const ccomplex_t IT_0040 = e_em*IT_0039;
    const ccomplex_t IT_0041 = IT_0000*IT_0010;
    const ccomplex_t IT_0042 = e_em*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0040 + (-3)*IT_0042);
    const ccomplex_t IT_0044 = (-0.166666666666667)*IT_0043;
    const ccomplex_t IT_0045 = IT_0036*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = (IT_0032 + -IT_0034)*(IT_0038 + -IT_0046);
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = pow(m_c, 2);
    const ccomplex_t IT_0051 = (-18)*IT_0050;
    const ccomplex_t IT_0052 = 18*s_34;
    const ccomplex_t IT_0053 = IT_0051 + IT_0052;
    const ccomplex_t IT_0054 = -(IT_0032 + -IT_0034)*(IT_0038 + -IT_0046);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0058 = 0.666666666666667*IT_0057;
    const ccomplex_t IT_0059 = (-1.33333333333333)*IT_0057;
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0065 = (-2)*IT_0064;
    const ccomplex_t IT_0066 = IT_0061*IT_0064;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = 0.5*IT_0063 + (-0.0833333333333333)*IT_0068;
    const ccomplex_t IT_0070 = cpow(IT_0069, 2);
    const ccomplex_t IT_0071 = s_34*IT_0016;
    const ccomplex_t IT_0072 = (-36)*IT_0071;
    const ccomplex_t IT_0073 = IT_0016*IT_0050;
    const ccomplex_t IT_0074 = (-36)*IT_0073;
    const ccomplex_t IT_0075 = s_13*s_14;
    const ccomplex_t IT_0076 = 72*IT_0075;
    const ccomplex_t IT_0077 = IT_0072 + IT_0074 + IT_0076;
    const ccomplex_t IT_0078 = 0.25*IT_0068;
    const ccomplex_t IT_0079 = IT_0033*IT_0038;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = IT_0033*IT_0046;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = (-6)*IT_0071;
    const ccomplex_t IT_0084 = (-6)*IT_0073;
    const ccomplex_t IT_0085 = 12*IT_0075;
    const ccomplex_t IT_0086 = IT_0083 + IT_0084 + IT_0085;
    const ccomplex_t IT_0087 = m_c*s_13;
    const ccomplex_t IT_0088 = (-18)*IT_0071;
    const ccomplex_t IT_0089 = 36*IT_0075;
    const ccomplex_t IT_0090 = IT_0088 + IT_0089;
    const ccomplex_t IT_0091 = 18*conj(IT_0049);
    const ccomplex_t IT_0092 = -s_14;
    const ccomplex_t IT_0093 = s_13 + IT_0092;
    const ccomplex_t IT_0094 = m_c*IT_0093;
    const ccomplex_t IT_0095 = 18*IT_0028;
    const ccomplex_t IT_0096 = IT_0080 + IT_0082;
    const ccomplex_t IT_0097 = 18*conj(IT_0028);
    const ccomplex_t IT_0098 = (-18)*IT_0049;
    const ccomplex_t IT_0099 = ((-12)*IT_0028 + 12*conj(IT_0028) + (-6)
      *IT_0049 + 6*conj(IT_0049) + (-6)*IT_0056 + 6*conj(IT_0056))*IT_0078 + 
      (conj(IT_0080) + conj(IT_0082))*IT_0095 + IT_0096*IT_0097 + IT_0069*((-36)
      *IT_0028 + 36*conj(IT_0028) + (-18)*IT_0056 + 18*conj(IT_0056) + IT_0091 +
       IT_0098);
    const ccomplex_t IT_0100 = m_c*s_14;
    const ccomplex_t IT_0101 = (-18)*IT_0082;
    const ccomplex_t IT_0102 = (-18)*IT_0073;
    const ccomplex_t IT_0103 = IT_0088 + IT_0089 + IT_0102;
    const ccomplex_t IT_0104 = (-18)*conj(IT_0049);
    const ccomplex_t IT_0105 = (-18)*conj(IT_0082);
    const ccomplex_t IT_0106 = -IT_0069;
    const ccomplex_t IT_0107 = IT_0028*conj(IT_0049)*IT_0053 + conj(IT_0028)*
      (IT_0028*(36*s_34 + (-36)*IT_0050) + IT_0053*(IT_0049 + IT_0056)) + 
      -IT_0070*IT_0077 + (-2)*IT_0078*(IT_0069*((-12)*IT_0071 + (-12)*IT_0073 +
       24*IT_0075) + 0.5*IT_0077*IT_0078 + 0.5*(IT_0080 + -conj(IT_0080) +
       IT_0082 + -conj(IT_0082))*IT_0086) + IT_0049*(conj(IT_0049)*IT_0052 + 18
      *conj(IT_0082)*IT_0087) + IT_0082*(conj(IT_0082)*IT_0090 + IT_0087*IT_0091
      ) + IT_0094*IT_0099 + conj(IT_0056)*(IT_0028*IT_0053 + IT_0052*IT_0056 +
       18*IT_0080*IT_0087 + IT_0050*IT_0098 + IT_0100*IT_0101) + conj(IT_0080)*
      (18*IT_0056*IT_0087 + IT_0080*IT_0090 + IT_0098*IT_0100 + IT_0073*IT_0101 
      + IT_0069*IT_0103) + IT_0080*(IT_0100*IT_0104 + IT_0073*IT_0105) + IT_0056
      *(IT_0050*IT_0104 + IT_0100*IT_0105) + IT_0103*(IT_0069*conj(IT_0082) +
       IT_0096*IT_0106);
    return create_ccomplex_return(IT_0107);
}

