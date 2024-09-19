#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_mu_anti_snu_mu_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_mu_anti_snu_mu_to_anti_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
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
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cos(alpha);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = cpow(IT_0001, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*m_W*e_em*(IT_0005*IT_0006 +
       IT_0007*IT_0008)*IT_0009;
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = IT_0009 + IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_W*e_em*(IT_0005*IT_0006 + 
      -IT_0007*IT_0008)*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + (-2)*IT_0018 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_W*e_em*(IT_0006*IT_0007 + 
      -IT_0005*IT_0008)*IT_0009;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_W*e_em*(IT_0006*IT_0007 +
       IT_0005*IT_0008)*IT_0014;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = s_12 + IT_0018;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0011;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = cpow(IT_0011, -1);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = IT_0033*IT_0036;
    const ccomplex_t IT_0038 = IT_0031*IT_0037;
    const ccomplex_t IT_0039 = (-2)*IT_0032;
    const ccomplex_t IT_0040 = (-0.5)*IT_0035;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = s_13 + s_23;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = IT_0038 + IT_0043;
    const ccomplex_t IT_0045 = IT_0036*IT_0039;
    const ccomplex_t IT_0046 = s_23*IT_0045;
    const ccomplex_t IT_0047 = IT_0033*IT_0040;
    const ccomplex_t IT_0048 = IT_0031*IT_0047;
    const ccomplex_t IT_0049 = -IT_0046 + (-2)*IT_0048;
    const ccomplex_t IT_0050 = IT_0044 + IT_0049;
    const ccomplex_t IT_0051 = IT_0030*IT_0050;
    const ccomplex_t IT_0052 = IT_0004 + -IT_0021 + -IT_0028 + IT_0051;
    const ccomplex_t IT_0053 = pow(m_W, -4);
    const ccomplex_t IT_0054 = pow(s_34, 2);
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = 2*IT_0032;
    const ccomplex_t IT_0057 = IT_0040*IT_0056;
    const ccomplex_t IT_0058 = IT_0036*IT_0056;
    const ccomplex_t IT_0059 = IT_0030*(IT_0047 + IT_0057 + -IT_0058);
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = s_13*s_24*s_34*IT_0053;
    const ccomplex_t IT_0062 = pow(m_W, -2);
    const ccomplex_t IT_0063 = s_14*s_24*IT_0062;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = s_13*s_23*IT_0062;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = s_12 + IT_0064 + IT_0066;
    const ccomplex_t IT_0068 = IT_0061 + IT_0067;
    const ccomplex_t IT_0069 = -IT_0030*(IT_0037 + -IT_0047 + -IT_0057);
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = s_14*s_23*s_34*IT_0053;
    const ccomplex_t IT_0072 = IT_0067 + IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0009;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = cpow(IT_0074, 2);
    const ccomplex_t IT_0076 = pow(m_W, 2);
    const ccomplex_t IT_0077 = cpow((-2)*s_23 + IT_0018 + IT_0076 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0075*IT_0078;
    const ccomplex_t IT_0080 = IT_0047 + IT_0057;
    const ccomplex_t IT_0081 = -IT_0037 + -IT_0058;
    const ccomplex_t IT_0082 = IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = IT_0030*IT_0082;
    const ccomplex_t IT_0084 = IT_0079 + -IT_0083;
    const ccomplex_t IT_0085 = s_23*s_24*s_34*IT_0053;
    const ccomplex_t IT_0086 = pow(s_24, 2);
    const ccomplex_t IT_0087 = IT_0062*IT_0086;
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = pow(s_23, 2);
    const ccomplex_t IT_0090 = IT_0062*IT_0089;
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = IT_0018 + IT_0085 + IT_0088 + IT_0091;
    const ccomplex_t IT_0093 = IT_0030*(IT_0047 + IT_0057);
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = pow(s_13, 2);
    const ccomplex_t IT_0096 = IT_0062*IT_0095;
    const ccomplex_t IT_0097 = pow(s_14, 2);
    const ccomplex_t IT_0098 = IT_0062*IT_0097;
    const ccomplex_t IT_0099 = -IT_0018;
    const ccomplex_t IT_0100 = s_13*s_14*s_34*IT_0053;
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = IT_0096 + IT_0098 + IT_0099 + IT_0101;
    const ccomplex_t IT_0103 = -IT_0052;
    const ccomplex_t IT_0104 = (s_12 + -IT_0063)*(s_12 + -IT_0065);
    const ccomplex_t IT_0105 = IT_0032*IT_0040;
    const ccomplex_t IT_0106 = IT_0032*IT_0036;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = IT_0105 + IT_0107;
    const ccomplex_t IT_0109 = IT_0030*IT_0108;
    const ccomplex_t IT_0110 = -IT_0079 + IT_0109;
    const ccomplex_t IT_0111 = s_14*s_34*IT_0062;
    const ccomplex_t IT_0112 = (s_12 + -IT_0065)*(s_13 + -IT_0111);
    const ccomplex_t IT_0113 = (IT_0096 + IT_0099)*(IT_0098 + IT_0099);
    const ccomplex_t IT_0114 = -conj(IT_0052);
    const ccomplex_t IT_0115 = -(s_12 + -IT_0063)*(IT_0018 + -IT_0096);
    const ccomplex_t IT_0116 = -(s_12 + -IT_0065)*(IT_0018 + -IT_0098);
    const ccomplex_t IT_0117 = -IT_0094;
    const ccomplex_t IT_0118 = -(IT_0018 + -IT_0087)*(IT_0018 + -IT_0096);
    const ccomplex_t IT_0119 = -IT_0118;
    const ccomplex_t IT_0120 = s_24*s_34*IT_0062;
    const ccomplex_t IT_0121 = (s_12 + -IT_0065)*(s_23 + -IT_0120);
    const ccomplex_t IT_0122 = -conj(IT_0094);
    const ccomplex_t IT_0123 = IT_0030*IT_0105;
    const ccomplex_t IT_0124 = -(IT_0018 + -IT_0090)*(IT_0018 + -IT_0098);
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = s_34*IT_0062;
    const ccomplex_t IT_0127 = s_13*s_34*IT_0062;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = s_14 + IT_0128;
    const ccomplex_t IT_0130 = IT_0126*IT_0129;
    const ccomplex_t IT_0131 = -(IT_0018 + -IT_0096)*(s_13 + -IT_0111);
    const ccomplex_t IT_0132 = -(IT_0018 + -IT_0096)*(s_23 + -IT_0120);
    const ccomplex_t IT_0133 = -IT_0060;
    const ccomplex_t IT_0134 = IT_0054*IT_0062;
    const ccomplex_t IT_0135 = -(IT_0018 + -IT_0096)*(IT_0076 + -IT_0134);
    const ccomplex_t IT_0136 = (s_12 + -IT_0065)*(IT_0076 + -IT_0134);
    const ccomplex_t IT_0137 = -conj(IT_0060);
    const ccomplex_t IT_0138 = (s_12 + -IT_0065)*(IT_0018 + -IT_0087);
    const ccomplex_t IT_0139 = -IT_0138;
    const ccomplex_t IT_0140 = (s_12 + -IT_0063)*(IT_0018 + -IT_0090);
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = -IT_0070;
    const ccomplex_t IT_0143 = s_23*s_34*IT_0062;
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = s_24 + IT_0144;
    const ccomplex_t IT_0146 = IT_0126*IT_0145;
    const ccomplex_t IT_0147 = (IT_0018 + -IT_0090)*(s_13 + -IT_0111);
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = (IT_0018 + -IT_0087)*(IT_0018 + -IT_0090);
    const ccomplex_t IT_0150 = -conj(IT_0070);
    const ccomplex_t IT_0151 = (IT_0018 + -IT_0090)*(s_23 + -IT_0120);
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = -conj(IT_0110);
    const ccomplex_t IT_0154 = (IT_0018 + -IT_0090)*(IT_0076 + -IT_0134);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = 2*IT_0052*(conj(IT_0052)*(1 + 0.5*IT_0055) +
       0.5*conj(IT_0060)*IT_0068 + 0.5*conj(IT_0070)*IT_0072 + 0.5*conj(IT_0084)
      *IT_0092) + conj(IT_0094)*IT_0102*IT_0103 + IT_0094*(conj(IT_0084)*IT_0104
       + conj(IT_0110)*IT_0112 + conj(IT_0094)*IT_0113 + IT_0102*IT_0114) + 
      (conj(IT_0060)*IT_0115 + conj(IT_0070)*IT_0116)*IT_0117 + IT_0060*(conj
      (IT_0052)*IT_0068 + conj(IT_0070)*IT_0104 + conj(IT_0060)*IT_0119 + conj
      (IT_0110)*IT_0121 + IT_0115*IT_0122) + IT_0070*(conj(IT_0052)*IT_0072 +
       conj(IT_0060)*IT_0104 + IT_0116*IT_0122 + IT_0112*conj(IT_0123) + conj
      (IT_0070)*IT_0125) + conj(IT_0123)*(IT_0103*IT_0130 + IT_0117*IT_0131 +
       IT_0132*IT_0133) + IT_0123*(conj(IT_0070)*IT_0112 + conj(IT_0084)*IT_0121
       + IT_0114*IT_0130 + IT_0122*IT_0131 + -conj(IT_0123)*IT_0135 + conj
      (IT_0110)*IT_0136 + IT_0132*IT_0137) + conj(IT_0084)*(IT_0133*IT_0139 +
       IT_0141*IT_0142) + conj(IT_0110)*(IT_0103*IT_0146 + IT_0142*IT_0148) +
       IT_0084*(conj(IT_0052)*IT_0092 + conj(IT_0094)*IT_0104 + IT_0121*conj
      (IT_0123) + IT_0137*IT_0139 + conj(IT_0084)*IT_0149 + IT_0141*IT_0150 +
       IT_0152*IT_0153) + IT_0110*(conj(IT_0094)*IT_0112 + conj(IT_0060)*IT_0121
       + conj(IT_0123)*IT_0136 + IT_0114*IT_0146 + IT_0148*IT_0150 + -conj
      (IT_0084)*IT_0152 + IT_0153*IT_0155);
    return create_ccomplex_return(IT_0156);
}

