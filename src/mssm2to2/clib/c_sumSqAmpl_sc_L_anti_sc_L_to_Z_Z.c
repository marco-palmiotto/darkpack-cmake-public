#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_sc_L_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_sc_L_to_Z_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -4);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = cpow(IT_0001, -2);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0003, -4);
    const ccomplex_t IT_0008 = cpow(IT_0001, 2);
    const ccomplex_t IT_0009 = tan(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, 2);
    const ccomplex_t IT_0011 = cpow(1 + IT_0010, -1);
    const ccomplex_t IT_0012 = (IT_0002 + (-0.666666666666667)*IT_0006 +
       0.111111111111111*IT_0007)*IT_0008*IT_0011;
    const ccomplex_t IT_0013 = IT_0000*IT_0012;
    const ccomplex_t IT_0014 = 9*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = 0.0555555555555556*IT_0015;
    const ccomplex_t IT_0017 = cos(beta);
    const ccomplex_t IT_0018 = cos(alpha);
    const ccomplex_t IT_0019 = IT_0001*IT_0007*IT_0017*IT_0018;
    const ccomplex_t IT_0020 = sin(beta);
    const ccomplex_t IT_0021 = sin(alpha);
    const ccomplex_t IT_0022 = IT_0001*IT_0007*IT_0020*IT_0021;
    const ccomplex_t IT_0023 = cpow(IT_0001, -3);
    const ccomplex_t IT_0024 = IT_0017*IT_0018*IT_0023;
    const ccomplex_t IT_0025 = IT_0020*IT_0021*IT_0023;
    const ccomplex_t IT_0026 = cpow(IT_0001, -1);
    const ccomplex_t IT_0027 = IT_0004*IT_0017*IT_0018*IT_0026;
    const ccomplex_t IT_0028 = IT_0004*IT_0020*IT_0021*IT_0026;
    const ccomplex_t IT_0029 = cpow(IT_0004 + IT_0005, -1);
    const ccomplex_t IT_0030 = (IT_0019 + IT_0022 + IT_0024 + IT_0025 + 2
      *IT_0027 + 2*IT_0028)*IT_0029;
    const ccomplex_t IT_0031 = m_W*e_em;
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = pow(m_W, -1);
    const ccomplex_t IT_0035 = pow(m_c, 2);
    const ccomplex_t IT_0036 = cpow(IT_0020, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*e_em*((-0.166666666666667)
      *m_W*(IT_0017*IT_0018 + -IT_0020*IT_0021)*(IT_0001*IT_0004 + (-3)*IT_0026)
       + IT_0021*IT_0026*IT_0034*IT_0035*IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = IT_0033*IT_0038;
    const ccomplex_t IT_0040 = pow(m_sc_L, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + (-2)*IT_0040 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0001*IT_0007*IT_0018*IT_0020;
    const ccomplex_t IT_0045 = IT_0001*IT_0007*IT_0017*IT_0021;
    const ccomplex_t IT_0046 = IT_0018*IT_0020*IT_0023;
    const ccomplex_t IT_0047 = IT_0017*IT_0021*IT_0023;
    const ccomplex_t IT_0048 = IT_0004*IT_0018*IT_0020*IT_0026;
    const ccomplex_t IT_0049 = IT_0004*IT_0017*IT_0021*IT_0026;
    const ccomplex_t IT_0050 = IT_0029*(IT_0044 + -IT_0045 + IT_0046 + 
      -IT_0047 + 2*IT_0048 + (-2)*IT_0049);
    const ccomplex_t IT_0051 = IT_0031*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*e_em*(0.166666666666667*m_W*
      (IT_0018*IT_0020 + IT_0017*IT_0021)*(IT_0001*IT_0004 + (-3)*IT_0026) +
       IT_0018*IT_0026*IT_0034*IT_0035*IT_0036);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0040 + -reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0016 + -IT_0043 + -IT_0058;
    const ccomplex_t IT_0060 = pow(m_Z, -4);
    const ccomplex_t IT_0061 = pow(s_34, 2);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = cpow(IT_0003, -1);
    const ccomplex_t IT_0064 = (0 + _Complex_I*-3)*e_em*(IT_0003*IT_0026 + -1.
      /3*IT_0001*IT_0063);
    const ccomplex_t IT_0065 = 0.166666666666667*IT_0064;
    const ccomplex_t IT_0066 = 0.333333333333333*IT_0064;
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = cpow(IT_0066, 2);
    const ccomplex_t IT_0069 = pow(m_Z, 2);
    const ccomplex_t IT_0070 = cpow((-2)*s_13 + IT_0069 + (0 + _Complex_I*1)
      *m_sc_L*Gamma_cl + reg_prop, -1);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = (IT_0067 + -IT_0068)*IT_0071;
    const ccomplex_t IT_0073 = pow(m_Z, -2);
    const ccomplex_t IT_0074 = pow(s_13, 2);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = pow(s_24, 2);
    const ccomplex_t IT_0077 = IT_0073*IT_0076;
    const ccomplex_t IT_0078 = -(IT_0040 + -IT_0075)*(IT_0040 + -IT_0077);
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = 3*IT_0072;
    const ccomplex_t IT_0081 = s_13*s_24*s_34*IT_0060;
    const ccomplex_t IT_0082 = s_14*s_24*IT_0073;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = s_13*s_23*IT_0073;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = s_12 + IT_0083 + IT_0085;
    const ccomplex_t IT_0087 = IT_0081 + IT_0086;
    const ccomplex_t IT_0088 = 3*IT_0059;
    const ccomplex_t IT_0089 = -(IT_0040 + -IT_0075)*(s_12 + -IT_0082);
    const ccomplex_t IT_0090 = IT_0067*IT_0071;
    const ccomplex_t IT_0091 = (-3)*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_23 + IT_0069 + (0 + _Complex_I*1)
      *m_sc_L*Gamma_cl + reg_prop, -1);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = IT_0067*IT_0093;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = (s_12 + -IT_0082)*(s_12 + -IT_0084);
    const ccomplex_t IT_0097 = s_14*s_34*IT_0073;
    const ccomplex_t IT_0098 = (s_12 + -IT_0084)*(s_13 + -IT_0097);
    const ccomplex_t IT_0099 = -IT_0090;
    const ccomplex_t IT_0100 = 3*IT_0099;
    const ccomplex_t IT_0101 = s_14*s_23*s_34*IT_0060;
    const ccomplex_t IT_0102 = IT_0086 + IT_0101;
    const ccomplex_t IT_0103 = pow(s_14, 2);
    const ccomplex_t IT_0104 = IT_0073*IT_0103;
    const ccomplex_t IT_0105 = -(s_12 + -IT_0084)*(IT_0040 + -IT_0104);
    const ccomplex_t IT_0106 = 3*conj(IT_0059);
    const ccomplex_t IT_0107 = (-3)*conj(IT_0090);
    const ccomplex_t IT_0108 = (IT_0067 + -IT_0068)*IT_0093;
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = 3*IT_0090;
    const ccomplex_t IT_0111 = s_24*s_34*IT_0073;
    const ccomplex_t IT_0112 = (s_12 + -IT_0084)*(s_23 + -IT_0111);
    const ccomplex_t IT_0113 = s_23*s_24*s_34*IT_0060;
    const ccomplex_t IT_0114 = -IT_0077;
    const ccomplex_t IT_0115 = pow(s_23, 2);
    const ccomplex_t IT_0116 = IT_0073*IT_0115;
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = IT_0040 + IT_0113 + IT_0114 + IT_0117;
    const ccomplex_t IT_0119 = (IT_0040 + -IT_0077)*(s_12 + -IT_0084);
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = (-3)*IT_0072;
    const ccomplex_t IT_0122 = (s_12 + -IT_0082)*(IT_0040 + -IT_0116);
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = (-3)*IT_0095;
    const ccomplex_t IT_0125 = (s_23 + -IT_0111)*(IT_0040 + -IT_0116);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = (-3)*IT_0108;
    const ccomplex_t IT_0128 = IT_0061*IT_0073;
    const ccomplex_t IT_0129 = (s_12 + -IT_0084)*(IT_0069 + -IT_0128);
    const ccomplex_t IT_0130 = (s_13 + -IT_0097)*(IT_0040 + -IT_0116);
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = (IT_0040 + -IT_0116)*(IT_0069 + -IT_0128);
    const ccomplex_t IT_0133 = -IT_0132;
    const ccomplex_t IT_0134 = s_34*IT_0073;
    const ccomplex_t IT_0135 = s_23*s_34*IT_0073;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = s_24 + IT_0136;
    const ccomplex_t IT_0138 = IT_0134*IT_0137;
    const ccomplex_t IT_0139 = (-3)*IT_0059;
    const ccomplex_t IT_0140 = -(IT_0040 + -IT_0075)*(s_13 + -IT_0097);
    const ccomplex_t IT_0141 = -(IT_0040 + -IT_0075)*(s_23 + -IT_0111);
    const ccomplex_t IT_0142 = s_13*s_34*IT_0073;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = s_14 + IT_0143;
    const ccomplex_t IT_0145 = IT_0134*IT_0144;
    const ccomplex_t IT_0146 = -IT_0040;
    const ccomplex_t IT_0147 = (IT_0075 + IT_0146)*(IT_0104 + IT_0146);
    const ccomplex_t IT_0148 = s_13*s_14*s_34*IT_0060;
    const ccomplex_t IT_0149 = IT_0075 + IT_0104 + IT_0146 + -IT_0148;
    const ccomplex_t IT_0150 = (-3)*conj(IT_0059);
    const ccomplex_t IT_0151 = 3*conj(IT_0090);
    const ccomplex_t IT_0152 = 3*conj(IT_0072);
    const ccomplex_t IT_0153 = 3*conj(IT_0099);
    const ccomplex_t IT_0154 = (-3)*conj(IT_0095);
    const ccomplex_t IT_0155 = -(IT_0040 + -IT_0075)*(IT_0069 + -IT_0128);
    const ccomplex_t IT_0156 = (-3)*conj(IT_0072);
    const ccomplex_t IT_0157 = (IT_0040 + -IT_0077)*(IT_0040 + -IT_0116);
    const ccomplex_t IT_0158 = 0.333333333333333*IT_0096;
    const ccomplex_t IT_0159 = 3*IT_0109*(0.333333333333333*IT_0106*IT_0118 + 
      -conj(IT_0108)*IT_0126 + 0.333333333333333*IT_0112*IT_0153 +
       0.333333333333333*IT_0123*IT_0154 + 0.333333333333333*IT_0120*IT_0156 +
       conj(IT_0109)*IT_0157 + IT_0151*IT_0158);
    const ccomplex_t IT_0160 = -(IT_0040 + -IT_0104)*(IT_0040 + -IT_0116);
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = 6*IT_0059*conj(IT_0059)*(1 + 0.5*IT_0062) +
       conj(IT_0072)*(IT_0079*IT_0080 + IT_0087*IT_0088 + IT_0089*IT_0091) +
       conj(IT_0095)*(IT_0080*IT_0096 + IT_0098*IT_0100 + IT_0088*IT_0102 +
       IT_0091*IT_0105) + IT_0072*(IT_0087*IT_0106 + IT_0089*IT_0107) + conj
      (IT_0109)*(IT_0096*IT_0110 + IT_0100*IT_0112 + IT_0088*IT_0118 + IT_0120
      *IT_0121 + IT_0123*IT_0124 + IT_0126*IT_0127) + conj(IT_0108)*(IT_0098
      *IT_0110 + IT_0080*IT_0112 + IT_0100*IT_0129 + IT_0124*IT_0131 + IT_0127
      *IT_0133 + IT_0138*IT_0139) + conj(IT_0099)*(IT_0091*IT_0140 + IT_0121
      *IT_0141 + IT_0139*IT_0145) + conj(IT_0090)*(IT_0110*IT_0147 + IT_0139
      *IT_0149) + IT_0090*IT_0149*IT_0150 + IT_0108*(IT_0138*IT_0150 + IT_0098
      *IT_0151 + IT_0112*IT_0152 + IT_0129*IT_0153 + IT_0131*IT_0154) + (-3)
      *IT_0099*((-0.333333333333333)*IT_0107*IT_0140 + (-0.333333333333333)
      *IT_0145*IT_0150 + conj(IT_0099)*IT_0155 + (-0.333333333333333)*IT_0141
      *IT_0156) + IT_0159 + 3*IT_0095*(0.333333333333333*IT_0102*IT_0106 +
       0.333333333333333*IT_0105*IT_0107 + 0.333333333333333*IT_0098*IT_0153 +
       IT_0152*IT_0158 + conj(IT_0095)*IT_0161);
    return create_ccomplex_return(IT_0162);
}

