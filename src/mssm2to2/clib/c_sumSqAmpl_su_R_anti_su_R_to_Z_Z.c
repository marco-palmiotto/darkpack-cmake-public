#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_su_R_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_su_R_to_Z_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
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
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -4);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = tan(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = cpow(1 + IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0000*IT_0002*IT_0004
      *IT_0007;
    const ccomplex_t IT_0009 = 0.888888888888889*IT_0008;
    const ccomplex_t IT_0010 = m_W*e_em;
    const ccomplex_t IT_0011 = cos(alpha);
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = IT_0002*IT_0003*IT_0011*IT_0012;
    const ccomplex_t IT_0014 = cos(beta);
    const ccomplex_t IT_0015 = sin(alpha);
    const ccomplex_t IT_0016 = IT_0002*IT_0003*IT_0014*IT_0015;
    const ccomplex_t IT_0017 = cpow(IT_0003, -3);
    const ccomplex_t IT_0018 = IT_0011*IT_0012*IT_0017;
    const ccomplex_t IT_0019 = IT_0014*IT_0015*IT_0017;
    const ccomplex_t IT_0020 = cpow(IT_0001, -2);
    const ccomplex_t IT_0021 = cpow(IT_0003, -1);
    const ccomplex_t IT_0022 = IT_0011*IT_0012*IT_0020*IT_0021;
    const ccomplex_t IT_0023 = IT_0014*IT_0015*IT_0020*IT_0021;
    const ccomplex_t IT_0024 = cpow(IT_0003, -2);
    const ccomplex_t IT_0025 = cpow(IT_0020 + IT_0024, -1);
    const ccomplex_t IT_0026 = (IT_0013 + -IT_0016 + IT_0018 + -IT_0019 + 2
      *IT_0022 + (-2)*IT_0023)*IT_0025;
    const ccomplex_t IT_0027 = IT_0010*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = pow(m_u, 2);
    const ccomplex_t IT_0031 = cpow(IT_0012, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*(IT_0011
      *IT_0012 + IT_0014*IT_0015)*IT_0020 + (-1.5)*IT_0011*IT_0021*IT_0029
      *IT_0030*IT_0031);
    const ccomplex_t IT_0033 = 0.666666666666667*IT_0032;
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = pow(m_su_R, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0035 + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0002*IT_0003*IT_0011*IT_0014;
    const ccomplex_t IT_0040 = IT_0002*IT_0003*IT_0012*IT_0015;
    const ccomplex_t IT_0041 = IT_0011*IT_0014*IT_0017;
    const ccomplex_t IT_0042 = IT_0012*IT_0015*IT_0017;
    const ccomplex_t IT_0043 = IT_0011*IT_0014*IT_0020*IT_0021;
    const ccomplex_t IT_0044 = IT_0012*IT_0015*IT_0020*IT_0021;
    const ccomplex_t IT_0045 = IT_0025*(IT_0039 + IT_0040 + IT_0041 + IT_0042 
      + 2*IT_0043 + 2*IT_0044);
    const ccomplex_t IT_0046 = IT_0010*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*(IT_0011
      *IT_0014 + -IT_0012*IT_0015)*IT_0020 + 1.5*IT_0015*IT_0021*IT_0029*IT_0030
      *IT_0031);
    const ccomplex_t IT_0049 = (-0.666666666666667)*IT_0048;
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + (-2)*IT_0035 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0009 + -IT_0038 + -IT_0053;
    const ccomplex_t IT_0055 = pow(m_Z, -4);
    const ccomplex_t IT_0056 = pow(s_34, 2);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = cpow(IT_0001, -1);
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0058;
    const ccomplex_t IT_0060 = 0.666666666666667*IT_0059;
    const ccomplex_t IT_0061 = 1.33333333333333*IT_0059;
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = cpow(IT_0061, 2);
    const ccomplex_t IT_0064 = pow(m_Z, 2);
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0064 + (0 + _Complex_I*1)
      *m_su_R*Gamma_ur + reg_prop, -1);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = (IT_0062 + -IT_0063)*IT_0066;
    const ccomplex_t IT_0068 = pow(m_Z, -2);
    const ccomplex_t IT_0069 = pow(s_13, 2);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = pow(s_24, 2);
    const ccomplex_t IT_0072 = IT_0068*IT_0071;
    const ccomplex_t IT_0073 = -(IT_0035 + -IT_0070)*(IT_0035 + -IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = 3*IT_0067;
    const ccomplex_t IT_0076 = s_13*s_24*s_34*IT_0055;
    const ccomplex_t IT_0077 = s_14*s_24*IT_0068;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = s_13*s_23*IT_0068;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = s_12 + IT_0078 + IT_0080;
    const ccomplex_t IT_0082 = IT_0076 + IT_0081;
    const ccomplex_t IT_0083 = 3*IT_0054;
    const ccomplex_t IT_0084 = -(IT_0035 + -IT_0070)*(s_12 + -IT_0077);
    const ccomplex_t IT_0085 = IT_0062*IT_0066;
    const ccomplex_t IT_0086 = (-3)*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_23 + IT_0064 + (0 + _Complex_I*1)
      *m_su_R*Gamma_ur + reg_prop, -1);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0062*IT_0088;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = (s_12 + -IT_0077)*(s_12 + -IT_0079);
    const ccomplex_t IT_0092 = s_14*s_34*IT_0068;
    const ccomplex_t IT_0093 = (s_12 + -IT_0079)*(s_13 + -IT_0092);
    const ccomplex_t IT_0094 = -IT_0085;
    const ccomplex_t IT_0095 = 3*IT_0094;
    const ccomplex_t IT_0096 = s_14*s_23*s_34*IT_0055;
    const ccomplex_t IT_0097 = IT_0081 + IT_0096;
    const ccomplex_t IT_0098 = pow(s_14, 2);
    const ccomplex_t IT_0099 = IT_0068*IT_0098;
    const ccomplex_t IT_0100 = -(s_12 + -IT_0079)*(IT_0035 + -IT_0099);
    const ccomplex_t IT_0101 = 3*conj(IT_0054);
    const ccomplex_t IT_0102 = (-3)*conj(IT_0085);
    const ccomplex_t IT_0103 = (IT_0062 + -IT_0063)*IT_0088;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = 3*IT_0085;
    const ccomplex_t IT_0106 = s_24*s_34*IT_0068;
    const ccomplex_t IT_0107 = (s_12 + -IT_0079)*(s_23 + -IT_0106);
    const ccomplex_t IT_0108 = s_23*s_24*s_34*IT_0055;
    const ccomplex_t IT_0109 = -IT_0072;
    const ccomplex_t IT_0110 = pow(s_23, 2);
    const ccomplex_t IT_0111 = IT_0068*IT_0110;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0035 + IT_0108 + IT_0109 + IT_0112;
    const ccomplex_t IT_0114 = (IT_0035 + -IT_0072)*(s_12 + -IT_0079);
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = (-3)*IT_0067;
    const ccomplex_t IT_0117 = (s_12 + -IT_0077)*(IT_0035 + -IT_0111);
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = (-3)*IT_0090;
    const ccomplex_t IT_0120 = (s_23 + -IT_0106)*(IT_0035 + -IT_0111);
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = (-3)*IT_0103;
    const ccomplex_t IT_0123 = IT_0056*IT_0068;
    const ccomplex_t IT_0124 = (s_12 + -IT_0079)*(IT_0064 + -IT_0123);
    const ccomplex_t IT_0125 = (s_13 + -IT_0092)*(IT_0035 + -IT_0111);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = (IT_0035 + -IT_0111)*(IT_0064 + -IT_0123);
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = s_34*IT_0068;
    const ccomplex_t IT_0130 = s_23*s_34*IT_0068;
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = s_24 + IT_0131;
    const ccomplex_t IT_0133 = IT_0129*IT_0132;
    const ccomplex_t IT_0134 = (-3)*IT_0054;
    const ccomplex_t IT_0135 = -(IT_0035 + -IT_0070)*(s_13 + -IT_0092);
    const ccomplex_t IT_0136 = -(IT_0035 + -IT_0070)*(s_23 + -IT_0106);
    const ccomplex_t IT_0137 = s_13*s_34*IT_0068;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = s_14 + IT_0138;
    const ccomplex_t IT_0140 = IT_0129*IT_0139;
    const ccomplex_t IT_0141 = -IT_0035;
    const ccomplex_t IT_0142 = (IT_0070 + IT_0141)*(IT_0099 + IT_0141);
    const ccomplex_t IT_0143 = s_13*s_14*s_34*IT_0055;
    const ccomplex_t IT_0144 = IT_0070 + IT_0099 + IT_0141 + -IT_0143;
    const ccomplex_t IT_0145 = (-3)*conj(IT_0054);
    const ccomplex_t IT_0146 = 3*conj(IT_0085);
    const ccomplex_t IT_0147 = 3*conj(IT_0067);
    const ccomplex_t IT_0148 = 3*conj(IT_0094);
    const ccomplex_t IT_0149 = (-3)*conj(IT_0090);
    const ccomplex_t IT_0150 = -(IT_0035 + -IT_0070)*(IT_0064 + -IT_0123);
    const ccomplex_t IT_0151 = (-3)*conj(IT_0067);
    const ccomplex_t IT_0152 = (IT_0035 + -IT_0072)*(IT_0035 + -IT_0111);
    const ccomplex_t IT_0153 = 0.333333333333333*IT_0091;
    const ccomplex_t IT_0154 = 3*IT_0104*(0.333333333333333*IT_0101*IT_0113 + 
      -conj(IT_0103)*IT_0121 + 0.333333333333333*IT_0107*IT_0148 +
       0.333333333333333*IT_0118*IT_0149 + 0.333333333333333*IT_0115*IT_0151 +
       conj(IT_0104)*IT_0152 + IT_0146*IT_0153);
    const ccomplex_t IT_0155 = -(IT_0035 + -IT_0099)*(IT_0035 + -IT_0111);
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = 6*IT_0054*conj(IT_0054)*(1 + 0.5*IT_0057) +
       conj(IT_0067)*(IT_0074*IT_0075 + IT_0082*IT_0083 + IT_0084*IT_0086) +
       conj(IT_0090)*(IT_0075*IT_0091 + IT_0093*IT_0095 + IT_0083*IT_0097 +
       IT_0086*IT_0100) + IT_0067*(IT_0082*IT_0101 + IT_0084*IT_0102) + conj
      (IT_0104)*(IT_0091*IT_0105 + IT_0095*IT_0107 + IT_0083*IT_0113 + IT_0115
      *IT_0116 + IT_0118*IT_0119 + IT_0121*IT_0122) + conj(IT_0103)*(IT_0093
      *IT_0105 + IT_0075*IT_0107 + IT_0095*IT_0124 + IT_0119*IT_0126 + IT_0122
      *IT_0128 + IT_0133*IT_0134) + conj(IT_0094)*(IT_0086*IT_0135 + IT_0116
      *IT_0136 + IT_0134*IT_0140) + conj(IT_0085)*(IT_0105*IT_0142 + IT_0134
      *IT_0144) + IT_0085*IT_0144*IT_0145 + IT_0103*(IT_0133*IT_0145 + IT_0093
      *IT_0146 + IT_0107*IT_0147 + IT_0124*IT_0148 + IT_0126*IT_0149) + (-3)
      *IT_0094*((-0.333333333333333)*IT_0102*IT_0135 + (-0.333333333333333)
      *IT_0140*IT_0145 + conj(IT_0094)*IT_0150 + (-0.333333333333333)*IT_0136
      *IT_0151) + IT_0154 + 3*IT_0090*(0.333333333333333*IT_0097*IT_0101 +
       0.333333333333333*IT_0100*IT_0102 + 0.333333333333333*IT_0093*IT_0148 +
       IT_0147*IT_0153 + conj(IT_0090)*IT_0156);
    return create_ccomplex_return(IT_0157);
}

