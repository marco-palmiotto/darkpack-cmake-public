#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_se_R_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_se_R_to_anti_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
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
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_W*e_em*(IT_0000*IT_0001 +
       IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -2);
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = pow(m_e, 2);
    const ccomplex_t IT_0011 = cpow(IT_0000, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*e_em*(m_W*(IT_0000*IT_0001 +
       -IT_0002*IT_0003)*IT_0004*IT_0008 + -IT_0001*IT_0005*IT_0009*IT_0010
      *IT_0011);
    const ccomplex_t IT_0013 = IT_0006*IT_0012;
    const ccomplex_t IT_0014 = pow(m_se_R, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + (-2)*IT_0014 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = cpow(s_12 + IT_0014 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = s_12 + IT_0014;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = 2*IT_0020;
    const ccomplex_t IT_0026 = IT_0021*IT_0025;
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = (-2)*IT_0020;
    const ccomplex_t IT_0029 = IT_0020*IT_0028;
    const ccomplex_t IT_0030 = s_13 + s_23;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = IT_0025*IT_0028;
    const ccomplex_t IT_0033 = s_23*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0027 + (-0.5)*IT_0031 + 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0024 + IT_0034;
    const ccomplex_t IT_0036 = IT_0019*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_W*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003)*IT_0005;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*e_em*(m_W*(IT_0001*IT_0002 +
       IT_0000*IT_0003)*IT_0004*IT_0008 + -IT_0003*IT_0005*IT_0009*IT_0010
      *IT_0011);
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*IT_0005*IT_0007;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0007, -1);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0048;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = IT_0047*IT_0050;
    const ccomplex_t IT_0052 = IT_0023*IT_0051;
    const ccomplex_t IT_0053 = (-2)*IT_0046;
    const ccomplex_t IT_0054 = (-2)*IT_0049;
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = s_23*IT_0055;
    const ccomplex_t IT_0057 = (-2)*IT_0052 + -IT_0056;
    const ccomplex_t IT_0058 = IT_0050*IT_0053;
    const ccomplex_t IT_0059 = IT_0030*IT_0058;
    const ccomplex_t IT_0060 = IT_0047*IT_0054;
    const ccomplex_t IT_0061 = IT_0023*IT_0060;
    const ccomplex_t IT_0062 = IT_0059 + IT_0061;
    const ccomplex_t IT_0063 = IT_0057 + IT_0062;
    const ccomplex_t IT_0064 = IT_0045*IT_0063;
    const ccomplex_t IT_0065 = -IT_0017 + IT_0036 + -IT_0043 + IT_0064;
    const ccomplex_t IT_0066 = pow(m_W, -4);
    const ccomplex_t IT_0067 = pow(s_34, 2);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = 2*IT_0046;
    const ccomplex_t IT_0070 = IT_0050*IT_0069;
    const ccomplex_t IT_0071 = IT_0051 + IT_0070;
    const ccomplex_t IT_0072 = -IT_0060;
    const ccomplex_t IT_0073 = IT_0071 + IT_0072;
    const ccomplex_t IT_0074 = IT_0045*IT_0073;
    const ccomplex_t IT_0075 = IT_0020*IT_0025;
    const ccomplex_t IT_0076 = IT_0022 + IT_0075;
    const ccomplex_t IT_0077 = -IT_0026;
    const ccomplex_t IT_0078 = IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = IT_0019*IT_0078;
    const ccomplex_t IT_0080 = -IT_0074 + 0.5*IT_0079;
    const ccomplex_t IT_0081 = s_14*s_23*s_34*IT_0066;
    const ccomplex_t IT_0082 = pow(m_W, -2);
    const ccomplex_t IT_0083 = s_14*s_24*IT_0082;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = s_13*s_23*IT_0082;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = s_12 + IT_0084 + IT_0086;
    const ccomplex_t IT_0088 = IT_0081 + IT_0087;
    const ccomplex_t IT_0089 = cpow(IT_0025, 2);
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = IT_0076 + IT_0090;
    const ccomplex_t IT_0092 = IT_0019*IT_0091;
    const ccomplex_t IT_0093 = IT_0054*IT_0069;
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = IT_0071 + IT_0094;
    const ccomplex_t IT_0096 = IT_0045*IT_0095;
    const ccomplex_t IT_0097 = 0.5*IT_0092 + -IT_0096;
    const ccomplex_t IT_0098 = s_13*s_24*s_34*IT_0066;
    const ccomplex_t IT_0099 = IT_0087 + IT_0098;
    const ccomplex_t IT_0100 = -IT_0026 + -IT_0089;
    const ccomplex_t IT_0101 = IT_0076 + IT_0100;
    const ccomplex_t IT_0102 = IT_0019*IT_0101;
    const ccomplex_t IT_0103 = -IT_0060 + -IT_0093;
    const ccomplex_t IT_0104 = IT_0071 + IT_0103;
    const ccomplex_t IT_0105 = IT_0045*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0102 + -IT_0105;
    const ccomplex_t IT_0107 = s_23*s_24*s_34*IT_0066;
    const ccomplex_t IT_0108 = pow(s_24, 2);
    const ccomplex_t IT_0109 = IT_0082*IT_0108;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = pow(s_23, 2);
    const ccomplex_t IT_0112 = IT_0082*IT_0111;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = IT_0014 + IT_0107 + IT_0110 + IT_0113;
    const ccomplex_t IT_0115 = IT_0019*IT_0076;
    const ccomplex_t IT_0116 = IT_0045*IT_0071;
    const ccomplex_t IT_0117 = 0.5*IT_0115 + -IT_0116;
    const ccomplex_t IT_0118 = pow(s_13, 2);
    const ccomplex_t IT_0119 = IT_0082*IT_0118;
    const ccomplex_t IT_0120 = pow(s_14, 2);
    const ccomplex_t IT_0121 = IT_0082*IT_0120;
    const ccomplex_t IT_0122 = -IT_0014;
    const ccomplex_t IT_0123 = s_13*s_14*s_34*IT_0066;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = IT_0119 + IT_0121 + IT_0122 + IT_0124;
    const ccomplex_t IT_0126 = -IT_0065;
    const ccomplex_t IT_0127 = (s_12 + -IT_0083)*(s_12 + -IT_0085);
    const ccomplex_t IT_0128 = IT_0046*IT_0050;
    const ccomplex_t IT_0129 = IT_0046*IT_0054;
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = IT_0128 + IT_0130;
    const ccomplex_t IT_0132 = IT_0045*IT_0131;
    const ccomplex_t IT_0133 = cpow(IT_0020, 2);
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = IT_0075 + IT_0134;
    const ccomplex_t IT_0136 = IT_0019*IT_0135;
    const ccomplex_t IT_0137 = IT_0132 + 0.5*IT_0136;
    const ccomplex_t IT_0138 = s_14*s_34*IT_0082;
    const ccomplex_t IT_0139 = (s_12 + -IT_0085)*(s_13 + -IT_0138);
    const ccomplex_t IT_0140 = (IT_0119 + IT_0122)*(IT_0121 + IT_0122);
    const ccomplex_t IT_0141 = -conj(IT_0065);
    const ccomplex_t IT_0142 = -(s_12 + -IT_0085)*(IT_0014 + -IT_0121);
    const ccomplex_t IT_0143 = -(s_12 + -IT_0083)*(IT_0014 + -IT_0119);
    const ccomplex_t IT_0144 = -IT_0117;
    const ccomplex_t IT_0145 = IT_0019*IT_0133;
    const ccomplex_t IT_0146 = IT_0045*IT_0128;
    const ccomplex_t IT_0147 = (-0.5)*IT_0145 + IT_0146;
    const ccomplex_t IT_0148 = -(IT_0014 + -IT_0112)*(IT_0014 + -IT_0121);
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = -conj(IT_0117);
    const ccomplex_t IT_0151 = s_24*s_34*IT_0082;
    const ccomplex_t IT_0152 = (s_12 + -IT_0085)*(s_23 + -IT_0151);
    const ccomplex_t IT_0153 = pow(m_W, 2);
    const ccomplex_t IT_0154 = IT_0067*IT_0082;
    const ccomplex_t IT_0155 = (s_12 + -IT_0085)*(IT_0153 + -IT_0154);
    const ccomplex_t IT_0156 = s_34*IT_0082;
    const ccomplex_t IT_0157 = s_13*s_34*IT_0082;
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = s_14 + IT_0158;
    const ccomplex_t IT_0160 = IT_0156*IT_0159;
    const ccomplex_t IT_0161 = -(IT_0014 + -IT_0119)*(s_13 + -IT_0138);
    const ccomplex_t IT_0162 = (s_12 + -IT_0083)*(IT_0014 + -IT_0112);
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = -IT_0080;
    const ccomplex_t IT_0165 = s_23*s_34*IT_0082;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = s_24 + IT_0166;
    const ccomplex_t IT_0168 = IT_0156*IT_0167;
    const ccomplex_t IT_0169 = (IT_0014 + -IT_0112)*(s_13 + -IT_0138);
    const ccomplex_t IT_0170 = -IT_0169;
    const ccomplex_t IT_0171 = -(IT_0014 + -IT_0119)*(s_23 + -IT_0151);
    const ccomplex_t IT_0172 = -IT_0147;
    const ccomplex_t IT_0173 = -(IT_0014 + -IT_0119)*(IT_0153 + -IT_0154);
    const ccomplex_t IT_0174 = -(IT_0014 + -IT_0109)*(IT_0014 + -IT_0119);
    const ccomplex_t IT_0175 = -IT_0174;
    const ccomplex_t IT_0176 = (s_12 + -IT_0085)*(IT_0014 + -IT_0109);
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = -conj(IT_0106);
    const ccomplex_t IT_0179 = (IT_0014 + -IT_0109)*(IT_0014 + -IT_0112);
    const ccomplex_t IT_0180 = -conj(IT_0080);
    const ccomplex_t IT_0181 = (IT_0014 + -IT_0112)*(s_23 + -IT_0151);
    const ccomplex_t IT_0182 = -IT_0181;
    const ccomplex_t IT_0183 = -conj(IT_0137);
    const ccomplex_t IT_0184 = (IT_0014 + -IT_0112)*(IT_0153 + -IT_0154);
    const ccomplex_t IT_0185 = -IT_0184;
    const ccomplex_t IT_0186 = 2*IT_0065*(conj(IT_0065)*(1 + 0.5*IT_0068) +
       0.5*conj(IT_0080)*IT_0088 + 0.5*conj(IT_0097)*IT_0099 + 0.5*conj(IT_0106)
      *IT_0114) + conj(IT_0117)*IT_0125*IT_0126 + IT_0117*(conj(IT_0106)*IT_0127
       + conj(IT_0137)*IT_0139 + conj(IT_0117)*IT_0140 + IT_0125*IT_0141) + 
      (conj(IT_0080)*IT_0142 + conj(IT_0097)*IT_0143)*IT_0144 + IT_0080*(conj
      (IT_0065)*IT_0088 + conj(IT_0097)*IT_0127 + IT_0139*conj(IT_0147) + conj
      (IT_0080)*IT_0149 + IT_0142*IT_0150) + IT_0147*(conj(IT_0080)*IT_0139 +
       conj(IT_0106)*IT_0152 + conj(IT_0137)*IT_0155 + IT_0141*IT_0160 + IT_0150
      *IT_0161) + conj(IT_0106)*IT_0163*IT_0164 + conj(IT_0137)*(IT_0126*IT_0168
       + IT_0164*IT_0170) + conj(IT_0097)*IT_0171*IT_0172 + conj(IT_0147)*
      (IT_0126*IT_0160 + IT_0144*IT_0161 + IT_0172*IT_0173) + IT_0097*(conj
      (IT_0065)*IT_0099 + conj(IT_0080)*IT_0127 + IT_0143*IT_0150 + conj(IT_0137
      )*IT_0152 + -conj(IT_0147)*IT_0171 + conj(IT_0097)*IT_0175 + IT_0177
      *IT_0178) + IT_0106*(conj(IT_0065)*IT_0114 + conj(IT_0117)*IT_0127 + conj
      (IT_0147)*IT_0152 + -conj(IT_0097)*IT_0177 + conj(IT_0106)*IT_0179 +
       IT_0163*IT_0180 + IT_0182*IT_0183) + IT_0137*(conj(IT_0117)*IT_0139 +
       conj(IT_0097)*IT_0152 + conj(IT_0147)*IT_0155 + IT_0141*IT_0168 + IT_0170
      *IT_0180 + IT_0178*IT_0182 + IT_0183*IT_0185);
    return create_ccomplex_return(IT_0186);
}

