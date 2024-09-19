#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_st_1_to_anti_b_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_st_1_to_anti_b_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*IT_0003;
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0007 = (-1.33333333333333)*IT_0006;
    const ccomplex_t IT_0008 = IT_0005*IT_0007;
    const ccomplex_t IT_0009 = (-0.666666666666667)*IT_0006;
    const ccomplex_t IT_0010 = IT_0005*IT_0009;
    const ccomplex_t IT_0011 = sin(beta);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = m_t*V_tb*V_u1*e_em*IT_0014*conj(U_st_10);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = IT_0009*IT_0019;
    const ccomplex_t IT_0022 = pow(m_b, 2);
    const ccomplex_t IT_0023 = pow(m_C_1, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (IT_0008 + -IT_0010 + -IT_0020 + IT_0021)
      *IT_0025;
    const ccomplex_t IT_0027 = pow(m_st_1, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0022 + IT_0027 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0006*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = -(IT_0005 + -IT_0019)*IT_0030;
    const ccomplex_t IT_0032 = 2*IT_0031;
    const ccomplex_t IT_0033 = 6*s_13;
    const ccomplex_t IT_0034 = s_14*IT_0033;
    const ccomplex_t IT_0035 = cos(beta);
    const ccomplex_t IT_0036 = cpow(IT_0035, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0001*IT_0014*IT_0036*conj(U_st_00);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0006*IT_0038;
    const ccomplex_t IT_0040 = IT_0028*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = 2*IT_0041;
    const ccomplex_t IT_0043 = m_b*m_C_1;
    const ccomplex_t IT_0044 = 6*s_14*IT_0043;
    const ccomplex_t IT_0045 = (IT_0008 + -IT_0020)*IT_0025;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = (-6)*s_13;
    const ccomplex_t IT_0048 = s_14*IT_0047;
    const ccomplex_t IT_0049 = (-0.333333333333333)*IT_0006;
    const ccomplex_t IT_0050 = cpow((-2)*s_12 + IT_0022 + -IT_0023 + -IT_0027 
      + -reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = -IT_0005;
    const ccomplex_t IT_0054 = IT_0019 + IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = -IT_0020;
    const ccomplex_t IT_0057 = IT_0008 + IT_0056;
    const ccomplex_t IT_0058 = IT_0025*IT_0057;
    const ccomplex_t IT_0059 = (-2)*IT_0055 + IT_0058;
    const ccomplex_t IT_0060 = s_34*IT_0047;
    const ccomplex_t IT_0061 = IT_0024*IT_0038;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0007*IT_0062;
    const ccomplex_t IT_0064 = (-6)*s_14*IT_0043;
    const ccomplex_t IT_0065 = IT_0038*IT_0049;
    const ccomplex_t IT_0066 = IT_0050*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0063 + (-2)*IT_0067;
    const ccomplex_t IT_0069 = (-6)*s_34*IT_0043;
    const ccomplex_t IT_0070 = -(IT_0007 + -IT_0009)*IT_0062;
    const ccomplex_t IT_0071 = IT_0023*IT_0033;
    const ccomplex_t IT_0072 = 6*IT_0023*IT_0043;
    const ccomplex_t IT_0073 = -IT_0041;
    const ccomplex_t IT_0074 = s_34*IT_0023;
    const ccomplex_t IT_0075 = 6*IT_0074;
    const ccomplex_t IT_0076 = -IT_0031;
    const ccomplex_t IT_0077 = m_b*s_14*m_C_1;
    const ccomplex_t IT_0078 = 6*IT_0077;
    const ccomplex_t IT_0079 = -IT_0067;
    const ccomplex_t IT_0080 = s_13*s_14;
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = IT_0074 + IT_0081;
    const ccomplex_t IT_0083 = 6*IT_0082;
    const ccomplex_t IT_0084 = IT_0023*IT_0047;
    const ccomplex_t IT_0085 = IT_0043*IT_0047;
    const ccomplex_t IT_0086 = pow(s_13, 2);
    const ccomplex_t IT_0087 = (-6)*IT_0086;
    const ccomplex_t IT_0088 = (-6)*IT_0023*IT_0043;
    const ccomplex_t IT_0089 = -(IT_0005 + -IT_0019)*IT_0052;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = (-6)*IT_0077;
    const ccomplex_t IT_0092 = IT_0026*IT_0064 + IT_0048*IT_0070 + IT_0042
      *IT_0071 + IT_0032*IT_0072 + IT_0073*IT_0075 + IT_0076*IT_0078 + IT_0079
      *IT_0083 + IT_0063*IT_0084 + IT_0059*IT_0085 + IT_0068*IT_0087 + IT_0046
      *IT_0088 + IT_0090*IT_0091;
    const ccomplex_t IT_0093 = IT_0026*IT_0048 + IT_0064*IT_0070 + IT_0032
      *IT_0071 + IT_0042*IT_0072 + IT_0075*IT_0076 + IT_0073*IT_0078 + IT_0046
      *IT_0084 + IT_0068*IT_0085 + IT_0059*IT_0087 + IT_0063*IT_0088 + IT_0083
      *IT_0090 + IT_0079*IT_0091;
    const ccomplex_t IT_0094 = IT_0033*IT_0043;
    const ccomplex_t IT_0095 = 6*IT_0086;
    const ccomplex_t IT_0096 = (-6)*IT_0074;
    const ccomplex_t IT_0097 = (-6)*IT_0082;
    const ccomplex_t IT_0098 = IT_0026*IT_0034 + IT_0044*IT_0070 + IT_0046
      *IT_0071 + IT_0063*IT_0072 + IT_0078*IT_0079 + IT_0032*IT_0084 + IT_0042
      *IT_0088 + IT_0073*IT_0091 + IT_0068*IT_0094 + IT_0059*IT_0095 + IT_0076
      *IT_0096 + IT_0090*IT_0097;
    const ccomplex_t IT_0099 = IT_0026*IT_0044 + IT_0034*IT_0070 + IT_0063
      *IT_0071 + IT_0046*IT_0072 + IT_0042*IT_0084 + IT_0032*IT_0088 + IT_0078
      *IT_0090 + IT_0076*IT_0091 + IT_0059*IT_0094 + IT_0068*IT_0095 + IT_0073
      *IT_0096 + IT_0079*IT_0097;
    const ccomplex_t IT_0100 = m_b*s_34*m_C_1;
    const ccomplex_t IT_0101 = 6*IT_0100;
    const ccomplex_t IT_0102 = (-6)*IT_0022;
    const ccomplex_t IT_0103 = IT_0043*IT_0102;
    const ccomplex_t IT_0104 = s_14*IT_0022;
    const ccomplex_t IT_0105 = (-6)*IT_0104;
    const ccomplex_t IT_0106 = s_13*s_34;
    const ccomplex_t IT_0107 = (-2)*IT_0106;
    const ccomplex_t IT_0108 = IT_0104 + IT_0107;
    const ccomplex_t IT_0109 = (-6)*IT_0108;
    const ccomplex_t IT_0110 = (-6)*IT_0100;
    const ccomplex_t IT_0111 = IT_0026*IT_0069 + IT_0060*IT_0070 + IT_0046
      *IT_0085 + IT_0063*IT_0087 + IT_0032*IT_0094 + IT_0042*IT_0095 + IT_0076
      *IT_0101 + s_13*IT_0068*IT_0102 + IT_0059*IT_0103 + IT_0079*IT_0105 +
       IT_0073*IT_0109 + IT_0090*IT_0110;
    const ccomplex_t IT_0112 = IT_0026*IT_0060 + IT_0069*IT_0070 + IT_0063
      *IT_0085 + IT_0046*IT_0087 + IT_0042*IT_0094 + IT_0032*IT_0095 + IT_0073
      *IT_0101 + s_13*IT_0059*IT_0102 + IT_0068*IT_0103 + IT_0090*IT_0105 +
       IT_0076*IT_0109 + IT_0079*IT_0110;
    const ccomplex_t IT_0113 = s_14*s_34;
    const ccomplex_t IT_0114 = 24*IT_0113;
    const ccomplex_t IT_0115 = (IT_0076 + IT_0090)*IT_0114;
    const ccomplex_t IT_0116 = IT_0073*IT_0114;
    const ccomplex_t IT_0117 = IT_0079*IT_0114;
    const ccomplex_t IT_0118 = conj(IT_0026)*(IT_0032*IT_0034 + IT_0042
      *IT_0044 + IT_0046*IT_0048 + IT_0059*IT_0060 + IT_0063*IT_0064 + IT_0068
      *IT_0069) + (IT_0034*IT_0042 + IT_0032*IT_0044 + IT_0048*IT_0063 + IT_0046
      *IT_0064 + IT_0060*IT_0068 + IT_0059*IT_0069)*conj(IT_0070) + conj(IT_0063
      )*IT_0092 + conj(IT_0046)*IT_0093 + conj(IT_0032)*IT_0098 + conj(IT_0042)
      *IT_0099 + conj(IT_0068)*IT_0111 + conj(IT_0059)*IT_0112 + conj(IT_0076)*
      (IT_0046*IT_0075 + IT_0063*IT_0078 + IT_0042*IT_0091 + IT_0032*IT_0096 +
       IT_0068*IT_0101 + IT_0059*IT_0109 + IT_0115) + conj(IT_0090)*(IT_0042
      *IT_0078 + IT_0046*IT_0083 + IT_0063*IT_0091 + IT_0032*IT_0097 + IT_0059
      *IT_0105 + IT_0068*IT_0110 + IT_0115) + conj(IT_0073)*(IT_0063*IT_0075 +
       IT_0046*IT_0078 + IT_0032*IT_0091 + IT_0042*IT_0096 + IT_0059*IT_0101 +
       IT_0068*IT_0109 + IT_0116 + IT_0117) + conj(IT_0079)*(IT_0032*IT_0078 +
       IT_0063*IT_0083 + IT_0046*IT_0091 + IT_0042*IT_0097 + IT_0068*IT_0105 +
       IT_0059*IT_0110 + IT_0116 + IT_0117);
    return create_ccomplex_return(IT_0118);
}

