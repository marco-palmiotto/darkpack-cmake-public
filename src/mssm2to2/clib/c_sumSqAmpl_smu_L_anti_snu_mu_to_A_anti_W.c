#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_snu_mu_to_A_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_snu_mu_to_A_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = (-0.5)*IT_0003;
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = pow(m_mu, 2);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, 2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0002*
      (IT_0005*IT_0006 + -m_W*(IT_0008 + -IT_0010));
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*m_W*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = pow(m_W, 2);
    const ccomplex_t IT_0016 = pow(m_smu_L, 2);
    const ccomplex_t IT_0017 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + IT_0015 + -IT_0016 + -IT_0017 
      + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0016 + -IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0002;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = s_13 + s_23;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = -IT_0023;
    const ccomplex_t IT_0031 = IT_0026*IT_0030;
    const ccomplex_t IT_0032 = s_23*IT_0031;
    const ccomplex_t IT_0033 = IT_0025*IT_0030;
    const ccomplex_t IT_0034 = -s_23;
    const ccomplex_t IT_0035 = s_12 + IT_0017;
    const ccomplex_t IT_0036 = IT_0034 + IT_0035;
    const ccomplex_t IT_0037 = IT_0033*IT_0036;
    const ccomplex_t IT_0038 = IT_0024*IT_0025;
    const ccomplex_t IT_0039 = (-0.5)*s_13 + (-0.5)*s_23 + 0.5*IT_0016 + 0.5
      *IT_0017;
    const ccomplex_t IT_0040 = s_12 + IT_0039;
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = -IT_0032 + -IT_0037 + 2*IT_0041;
    const ccomplex_t IT_0043 = IT_0029 + IT_0042;
    const ccomplex_t IT_0044 = IT_0022*IT_0043;
    const ccomplex_t IT_0045 = IT_0004 + -IT_0020 + IT_0044;
    const ccomplex_t IT_0046 = pow(m_W, -2);
    const ccomplex_t IT_0047 = pow(s_34, 2);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = 2*IT_0025;
    const ccomplex_t IT_0050 = IT_0024*IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0038 + IT_0051;
    const ccomplex_t IT_0053 = IT_0022*IT_0052;
    const ccomplex_t IT_0054 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_smu_L
      *Gamma_ml + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0050*IT_0055;
    const ccomplex_t IT_0057 = -IT_0053 + 0.5*IT_0056;
    const ccomplex_t IT_0058 = IT_0030*IT_0049;
    const ccomplex_t IT_0059 = IT_0022*(IT_0033 + -IT_0038 + IT_0050 + 
      -IT_0058);
    const ccomplex_t IT_0060 = (IT_0038 + -IT_0050)*IT_0055;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = -s_13;
    const ccomplex_t IT_0063 = IT_0057*(IT_0016*conj(IT_0057) + s_12*conj
      (IT_0059)) + IT_0059*(s_12*conj(IT_0057) + IT_0017*conj(IT_0059)) + conj
      (IT_0061)*(IT_0045 + -IT_0034*IT_0059 + -IT_0057*IT_0062) + IT_0061*(conj
      (IT_0045) + -IT_0034*conj(IT_0059) + -conj(IT_0057)*IT_0062);
    const ccomplex_t IT_0064 = pow(s_14, 2);
    const ccomplex_t IT_0065 = IT_0046*IT_0064;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0016 + IT_0066;
    const ccomplex_t IT_0068 = (-2)*IT_0025;
    const ccomplex_t IT_0069 = IT_0024*IT_0068;
    const ccomplex_t IT_0070 = IT_0038 + IT_0069;
    const ccomplex_t IT_0071 = IT_0022*IT_0070;
    const ccomplex_t IT_0072 = (-0.5)*IT_0056 + IT_0071;
    const ccomplex_t IT_0073 = IT_0022*(IT_0033 + -IT_0038 + -IT_0069);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = s_12*conj(IT_0072);
    const ccomplex_t IT_0076 = IT_0017*conj(IT_0074);
    const ccomplex_t IT_0077 = IT_0075 + IT_0076;
    const ccomplex_t IT_0078 = -IT_0069;
    const ccomplex_t IT_0079 = IT_0027 + IT_0078;
    const ccomplex_t IT_0080 = IT_0022*IT_0079;
    const ccomplex_t IT_0081 = IT_0052*IT_0055;
    const ccomplex_t IT_0082 = IT_0080 + (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = s_13*IT_0072;
    const ccomplex_t IT_0084 = s_23*IT_0074;
    const ccomplex_t IT_0085 = IT_0083 + IT_0084;
    const ccomplex_t IT_0086 = s_13*conj(IT_0072);
    const ccomplex_t IT_0087 = s_23*conj(IT_0074);
    const ccomplex_t IT_0088 = IT_0086 + IT_0087;
    const ccomplex_t IT_0089 = s_12*IT_0072;
    const ccomplex_t IT_0090 = IT_0016*IT_0072;
    const ccomplex_t IT_0091 = s_14*s_24*IT_0046;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = s_12 + IT_0092;
    const ccomplex_t IT_0094 = IT_0030*IT_0068;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0070 + IT_0095;
    const ccomplex_t IT_0097 = IT_0022*IT_0096;
    const ccomplex_t IT_0098 = -IT_0058;
    const ccomplex_t IT_0099 = IT_0050 + IT_0098;
    const ccomplex_t IT_0100 = IT_0055*IT_0099;
    const ccomplex_t IT_0101 = IT_0097 + (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = s_12*conj(IT_0074);
    const ccomplex_t IT_0103 = s_13*conj(IT_0082);
    const ccomplex_t IT_0104 = IT_0016*conj(IT_0072);
    const ccomplex_t IT_0105 = conj(IT_0045) + IT_0102 + IT_0103 + IT_0104;
    const ccomplex_t IT_0106 = IT_0022*(IT_0033 + -IT_0038 + -IT_0069 +
       IT_0094);
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = s_23*IT_0082;
    const ccomplex_t IT_0109 = IT_0017*IT_0074 + IT_0089 + IT_0108;
    const ccomplex_t IT_0110 = s_23*conj(IT_0082);
    const ccomplex_t IT_0111 = IT_0077 + IT_0110;
    const ccomplex_t IT_0112 = IT_0027 + IT_0094;
    const ccomplex_t IT_0113 = -IT_0031 + -IT_0069;
    const ccomplex_t IT_0114 = IT_0112 + IT_0113;
    const ccomplex_t IT_0115 = IT_0022*IT_0114;
    const ccomplex_t IT_0116 = IT_0033 + IT_0050;
    const ccomplex_t IT_0117 = -IT_0038 + -IT_0058;
    const ccomplex_t IT_0118 = IT_0116 + IT_0117;
    const ccomplex_t IT_0119 = IT_0055*IT_0118;
    const ccomplex_t IT_0120 = IT_0115 + 0.5*IT_0119;
    const ccomplex_t IT_0121 = s_13*IT_0082;
    const ccomplex_t IT_0122 = IT_0093*(conj(IT_0045)*IT_0074 + IT_0045*conj
      (IT_0074) + IT_0101*IT_0105 + conj(IT_0107)*IT_0109 + IT_0107*IT_0111 +
       IT_0088*IT_0120 + IT_0085*conj(IT_0120) + conj(IT_0101)*(IT_0045 + s_12
      *IT_0074 + IT_0090 + IT_0121));
    const ccomplex_t IT_0123 = s_14*s_34*IT_0046;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = s_13 + IT_0124;
    const ccomplex_t IT_0126 = (conj(IT_0045)*IT_0082 + IT_0045*conj(IT_0082) 
      + conj(IT_0061)*IT_0085 + IT_0061*IT_0088 + IT_0057*IT_0105 + conj(IT_0059
      )*IT_0109 + IT_0059*IT_0111 + conj(IT_0057)*(IT_0045 + s_12*IT_0074 +
       IT_0090 + IT_0121))*IT_0125;
    const ccomplex_t IT_0127 = s_24*s_34*IT_0046;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = s_23 + IT_0128;
    const ccomplex_t IT_0130 = s_13*IT_0101;
    const ccomplex_t IT_0131 = s_23*IT_0107;
    const ccomplex_t IT_0132 = IT_0130 + IT_0131;
    const ccomplex_t IT_0133 = s_13*conj(IT_0101);
    const ccomplex_t IT_0134 = s_23*conj(IT_0107);
    const ccomplex_t IT_0135 = IT_0133 + IT_0134;
    const ccomplex_t IT_0136 = s_12*IT_0101;
    const ccomplex_t IT_0137 = s_12*conj(IT_0101);
    const ccomplex_t IT_0138 = IT_0016*IT_0101;
    const ccomplex_t IT_0139 = conj(IT_0045)*IT_0120 + IT_0045*conj(IT_0120) +
       IT_0057*(IT_0016*conj(IT_0101) + s_12*conj(IT_0107) + s_13*conj(IT_0120))
       + conj(IT_0061)*IT_0132 + IT_0061*IT_0135 + conj(IT_0059)*(IT_0045 +
       IT_0017*IT_0107 + s_23*IT_0120 + IT_0136) + IT_0059*(conj(IT_0045) +
       IT_0017*conj(IT_0107) + s_23*conj(IT_0120) + IT_0137) + conj(IT_0057)*
      (s_12*IT_0107 + s_13*IT_0120 + IT_0138);
    const ccomplex_t IT_0140 = pow(s_24, 2);
    const ccomplex_t IT_0141 = IT_0046*IT_0140;
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = IT_0017 + IT_0142;
    const ccomplex_t IT_0144 = 3*IT_0045*conj(IT_0045) + -IT_0048*IT_0063 +
       IT_0067*(conj(IT_0045)*IT_0072 + IT_0074*IT_0077 + conj(IT_0082)*IT_0085 
      + IT_0082*IT_0088 + conj(IT_0074)*IT_0089 + conj(IT_0072)*(IT_0045 +
       IT_0090)) + IT_0122 + IT_0126 + IT_0129*IT_0139 + (conj(IT_0120)*IT_0132 
      + IT_0120*IT_0135 + conj(IT_0107)*(IT_0045 + IT_0136) + IT_0107*(conj
      (IT_0045) + IT_0017*conj(IT_0107) + IT_0137) + conj(IT_0101)*IT_0138)
      *IT_0143;
    return create_ccomplex_return(IT_0144);
}

