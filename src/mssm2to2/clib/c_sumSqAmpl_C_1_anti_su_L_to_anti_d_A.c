#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_su_L_to_anti_d_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_su_L_to_anti_d_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = s_14*s_34;
    const ccomplex_t IT_0001 = pow(m_d, 2);
    const ccomplex_t IT_0002 = pow(m_C_1, 2);
    const ccomplex_t IT_0003 = pow(m_su_L, 2);
    const ccomplex_t IT_0004 = cpow((-2)*s_12 + IT_0001 + -IT_0002 + -IT_0003 
      + -reg_prop, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp1*IT_0006;
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0010 = (-0.333333333333333)*IT_0009;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = IT_0004*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = IT_0008*IT_0009;
    const ccomplex_t IT_0015 = cpow((-2)*s_23 + IT_0001 + IT_0003 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = (-1.33333333333333)*IT_0009;
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0008*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = s_34*IT_0002;
    const ccomplex_t IT_0025 = s_13*s_14;
    const ccomplex_t IT_0026 = (-2)*IT_0025;
    const ccomplex_t IT_0027 = IT_0024 + IT_0026;
    const ccomplex_t IT_0028 = 6*conj(IT_0013);
    const ccomplex_t IT_0029 = 6*conj(IT_0017);
    const ccomplex_t IT_0030 = (-0.666666666666667)*IT_0009;
    const ccomplex_t IT_0031 = pow(m_W, -1);
    const ccomplex_t IT_0032 = cos(beta);
    const ccomplex_t IT_0033 = cpow(IT_0032, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0006*IT_0031*IT_0033;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = IT_0019*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -(IT_0018 + -IT_0030)*IT_0037;
    const ccomplex_t IT_0039 = IT_0009*IT_0035;
    const ccomplex_t IT_0040 = IT_0015*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = 2*IT_0041;
    const ccomplex_t IT_0043 = 6*s_13;
    const ccomplex_t IT_0044 = s_14*IT_0043;
    const ccomplex_t IT_0045 = (-2)*IT_0017;
    const ccomplex_t IT_0046 = m_d*m_C_1;
    const ccomplex_t IT_0047 = 6*s_14*IT_0046;
    const ccomplex_t IT_0048 = IT_0018*IT_0037;
    const ccomplex_t IT_0049 = (-6)*s_13;
    const ccomplex_t IT_0050 = s_14*IT_0049;
    const ccomplex_t IT_0051 = IT_0010*IT_0035;
    const ccomplex_t IT_0052 = IT_0004*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = -IT_0048 + (-2)*IT_0053;
    const ccomplex_t IT_0055 = s_34*IT_0049;
    const ccomplex_t IT_0056 = (-6)*s_14*IT_0046;
    const ccomplex_t IT_0057 = 2*IT_0013 + IT_0022;
    const ccomplex_t IT_0058 = (-6)*s_34*IT_0046;
    const ccomplex_t IT_0059 = -IT_0021*(IT_0018 + -IT_0030);
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = m_d*s_14*m_C_1;
    const ccomplex_t IT_0062 = 6*IT_0017;
    const ccomplex_t IT_0063 = (-6)*IT_0013;
    const ccomplex_t IT_0064 = IT_0062 + IT_0063;
    const ccomplex_t IT_0065 = IT_0002*IT_0043;
    const ccomplex_t IT_0066 = 6*IT_0002*IT_0046;
    const ccomplex_t IT_0067 = -IT_0041;
    const ccomplex_t IT_0068 = 6*IT_0024;
    const ccomplex_t IT_0069 = -IT_0053;
    const ccomplex_t IT_0070 = 6*IT_0027;
    const ccomplex_t IT_0071 = IT_0002*IT_0049;
    const ccomplex_t IT_0072 = IT_0046*IT_0049;
    const ccomplex_t IT_0073 = pow(s_13, 2);
    const ccomplex_t IT_0074 = (-6)*IT_0073;
    const ccomplex_t IT_0075 = (-6)*IT_0002*IT_0046;
    const ccomplex_t IT_0076 = IT_0038*IT_0050 + IT_0056*IT_0060 + IT_0061
      *IT_0064 + IT_0042*IT_0065 + IT_0045*IT_0066 + IT_0067*IT_0068 + IT_0069
      *IT_0070 + IT_0048*IT_0071 + IT_0057*IT_0072 + IT_0054*IT_0074 + IT_0023
      *IT_0075;
    const ccomplex_t IT_0077 = 6*IT_0013;
    const ccomplex_t IT_0078 = 6*IT_0061;
    const ccomplex_t IT_0079 = (-6)*IT_0061;
    const ccomplex_t IT_0080 = IT_0038*IT_0056 + IT_0050*IT_0060 + IT_0024
      *IT_0062 + IT_0045*IT_0065 + IT_0042*IT_0066 + IT_0023*IT_0071 + IT_0054
      *IT_0072 + IT_0057*IT_0074 + IT_0048*IT_0075 + IT_0027*IT_0077 + IT_0067
      *IT_0078 + IT_0069*IT_0079;
    const ccomplex_t IT_0081 = m_d*s_34*m_C_1;
    const ccomplex_t IT_0082 = IT_0043*IT_0046;
    const ccomplex_t IT_0083 = 6*IT_0073;
    const ccomplex_t IT_0084 = (-6)*IT_0001;
    const ccomplex_t IT_0085 = IT_0046*IT_0084;
    const ccomplex_t IT_0086 = s_14*IT_0001;
    const ccomplex_t IT_0087 = (-6)*IT_0086;
    const ccomplex_t IT_0088 = s_13*s_34;
    const ccomplex_t IT_0089 = (-2)*IT_0088;
    const ccomplex_t IT_0090 = IT_0086 + IT_0089;
    const ccomplex_t IT_0091 = (-6)*IT_0090;
    const ccomplex_t IT_0092 = IT_0038*IT_0055 + IT_0058*IT_0060 + IT_0023
      *IT_0072 + IT_0048*IT_0074 + IT_0064*IT_0081 + IT_0045*IT_0082 + IT_0042
      *IT_0083 + s_13*IT_0054*IT_0084 + IT_0057*IT_0085 + IT_0069*IT_0087 +
       IT_0067*IT_0091;
    const ccomplex_t IT_0093 = (-6)*IT_0017;
    const ccomplex_t IT_0094 = IT_0038*IT_0047 + IT_0044*IT_0060 + IT_0027
      *IT_0063 + IT_0023*IT_0065 + IT_0048*IT_0066 + IT_0045*IT_0071 + IT_0042
      *IT_0075 + IT_0069*IT_0078 + IT_0067*IT_0079 + IT_0054*IT_0082 + IT_0057
      *IT_0083 + IT_0024*IT_0093;
    const ccomplex_t IT_0095 = 6*IT_0081;
    const ccomplex_t IT_0096 = (-6)*IT_0081;
    const ccomplex_t IT_0097 = IT_0038*IT_0058 + IT_0055*IT_0060 + IT_0048
      *IT_0072 + IT_0023*IT_0074 + IT_0042*IT_0082 + IT_0045*IT_0083 + s_13
      *IT_0057*IT_0084 + IT_0054*IT_0085 + IT_0063*IT_0086 + IT_0090*IT_0093 +
       IT_0067*IT_0095 + IT_0069*IT_0096;
    const ccomplex_t IT_0098 = (-6)*IT_0024;
    const ccomplex_t IT_0099 = (-6)*IT_0027;
    const ccomplex_t IT_0100 = IT_0038*IT_0044 + IT_0047*IT_0060 + IT_0048
      *IT_0065 + IT_0023*IT_0066 + IT_0042*IT_0071 + IT_0045*IT_0075 + IT_0057
      *IT_0082 + IT_0054*IT_0083 + IT_0061*(IT_0077 + IT_0093) + IT_0067*IT_0098
       + IT_0069*IT_0099;
    const ccomplex_t IT_0101 = (-6)*conj(IT_0017);
    const ccomplex_t IT_0102 = (-6)*conj(IT_0013);
    const ccomplex_t IT_0103 = IT_0054*IT_0081;
    const ccomplex_t IT_0104 = IT_0067 + IT_0069;
    const ccomplex_t IT_0105 = 24*IT_0000;
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = IT_0000*(IT_0013 + IT_0017)*(24*conj(IT_0013) +
       24*conj(IT_0017)) + IT_0023*(IT_0027*IT_0028 + IT_0024*IT_0029) + conj
      (IT_0038)*(IT_0042*IT_0044 + IT_0045*IT_0047 + IT_0048*IT_0050 + IT_0054
      *IT_0055 + IT_0023*IT_0056 + IT_0057*IT_0058) + (IT_0044*IT_0045 + IT_0042
      *IT_0047 + IT_0023*IT_0050 + IT_0048*IT_0056 + IT_0055*IT_0057 + IT_0054
      *IT_0058)*conj(IT_0060) + conj(IT_0048)*IT_0076 + conj(IT_0023)*IT_0080 +
       conj(IT_0054)*IT_0092 + conj(IT_0045)*IT_0094 + conj(IT_0057)*IT_0097 +
       conj(IT_0042)*IT_0100 + (IT_0024*IT_0045 + IT_0057*IT_0090)*IT_0101 +
       IT_0061*(IT_0042*(IT_0028 + IT_0101) + IT_0048*(IT_0029 + IT_0102)) +
       IT_0102*(IT_0027*IT_0045 + IT_0057*IT_0086 + IT_0103) + IT_0029*IT_0103 +
       conj(IT_0067)*(IT_0048*IT_0068 + IT_0023*IT_0078 + IT_0045*IT_0079 +
       IT_0054*IT_0091 + IT_0057*IT_0095 + IT_0042*IT_0098 + IT_0106) + conj
      (IT_0069)*(IT_0048*IT_0070 + IT_0045*IT_0078 + IT_0023*IT_0079 + IT_0054
      *IT_0087 + IT_0057*IT_0096 + IT_0042*IT_0099 + IT_0106);
    return create_ccomplex_return(IT_0107);
}

