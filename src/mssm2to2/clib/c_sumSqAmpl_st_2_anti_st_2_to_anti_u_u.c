#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_anti_u_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_anti_u_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_01*conj(U_st_01) + IT_0001*IT_0002*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0005 = 1.33333333333333*IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0002;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0003;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0007 + (-3)*IT_0009);
    const ccomplex_t IT_0011 = (-0.166666666666667)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_st_2, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0005*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0019 = m_u*IT_0018;
    const ccomplex_t IT_0020 = m_u*IT_0005;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0022 = (-0.666666666666667)*IT_0021;
    const ccomplex_t IT_0023 = IT_0013*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (IT_0019 + -IT_0020)*(IT_0015 + -IT_0024);
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = m_u*s_13;
    const ccomplex_t IT_0028 = pow(m_W, -1);
    const ccomplex_t IT_0029 = m_t*e_em*mu_h*IT_0028*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0030 = IT_0003*IT_0029;
    const ccomplex_t IT_0031 = m_t*e_em*mu_h*IT_0028*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0032 = IT_0003*IT_0031;
    const ccomplex_t IT_0033 = cos(beta);
    const ccomplex_t IT_0034 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0038 = IT_0033*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0030 + 0.5*IT_0032 + (-0.5)*IT_0036 +
       0.5*IT_0039;
    const ccomplex_t IT_0041 = sin(beta);
    const ccomplex_t IT_0042 = cpow(IT_0041, -1);
    const ccomplex_t IT_0043 = m_u*e_em*IT_0003*IT_0028*IT_0033*IT_0042;
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0040*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_12 + (-2)*IT_0012 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = m_u*e_em*IT_0003*IT_0028;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0003*IT_0033*IT_0042;
    const ccomplex_t IT_0052 = IT_0029*IT_0051;
    const ccomplex_t IT_0053 = IT_0031*IT_0051;
    const ccomplex_t IT_0054 = IT_0034*IT_0041;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = IT_0037*IT_0041;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0052 + (-0.5)*IT_0053 + (-0.5)*IT_0055 +
       0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0050*IT_0058;
    const ccomplex_t IT_0060 = pow(m_Z, 2);
    const ccomplex_t IT_0061 = cpow((-2)*s_12 + (-2)*IT_0012 + IT_0060 + 
      -reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = -IT_0048 + -IT_0063;
    const ccomplex_t IT_0065 = s_13 + s_14;
    const ccomplex_t IT_0066 = m_u*IT_0065;
    const ccomplex_t IT_0067 = s_34*IT_0012;
    const ccomplex_t IT_0068 = (-18)*IT_0067;
    const ccomplex_t IT_0069 = s_13*s_14;
    const ccomplex_t IT_0070 = 36*IT_0069;
    const ccomplex_t IT_0071 = IT_0068 + IT_0070;
    const ccomplex_t IT_0072 = cos(alpha);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_u*e_em*IT_0003*IT_0028
      *IT_0042*IT_0072;
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = m_t*mu_h*IT_0028;
    const ccomplex_t IT_0076 = sin(alpha);
    const ccomplex_t IT_0077 = IT_0033*IT_0076;
    const ccomplex_t IT_0078 = IT_0041*IT_0072;
    const ccomplex_t IT_0079 = IT_0077 + IT_0078;
    const ccomplex_t IT_0080 = cpow(IT_0000, -2);
    const ccomplex_t IT_0081 = pow(m_t, 2);
    const ccomplex_t IT_0082 = 2*U_st_01*(((-0.5)*m_W*IT_0079*(IT_0003 + -1./3
      *IT_0002*IT_0080) + IT_0003*IT_0028*IT_0042*IT_0072*IT_0081)*conj(U_st_01)
       + 1./2*IT_0003*IT_0042*IT_0075*IT_0076*conj(U_st_11)) + U_st_11*(IT_0003
      *IT_0042*IT_0075*IT_0076*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0079*IT_0080 + (-1.5)*IT_0003*IT_0028*IT_0042*IT_0072*IT_0081)*conj
      (U_st_11));
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(e_em*IT_0082 +
       1.4142135623731*A_t*IT_0072*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = IT_0074*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*m_u*e_em*IT_0003*IT_0028
      *IT_0042*IT_0076;
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = IT_0033*IT_0072;
    const ccomplex_t IT_0092 = -IT_0041*IT_0076 + IT_0091;
    const ccomplex_t IT_0093 = (-2)*U_st_01*((IT_0003*IT_0028*IT_0042*IT_0076
      *IT_0081 + 0.5*m_W*(IT_0003 + -1./3*IT_0002*IT_0080)*IT_0092)*conj(U_st_01
      ) + -1./2*IT_0003*IT_0042*IT_0072*IT_0075*conj(U_st_11)) + U_st_11*
      (IT_0003*IT_0042*IT_0072*IT_0075*conj(U_st_01) + (-1.33333333333333)*(1.5
      *IT_0003*IT_0028*IT_0042*IT_0076*IT_0081 + m_W*IT_0002*IT_0080*IT_0092)
      *conj(U_st_11));
    const ccomplex_t IT_0094 = (0 + _Complex_I*-1)*(e_em*IT_0093 + (
      -1.4142135623731)*A_t*IT_0076*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0090*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_12 + (-2)*IT_0012 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = -IT_0088 + -IT_0099;
    const ccomplex_t IT_0101 = pow(m_u, 2);
    const ccomplex_t IT_0102 = 18*s_34;
    const ccomplex_t IT_0103 = (-18)*IT_0101;
    const ccomplex_t IT_0104 = IT_0102 + IT_0103;
    const ccomplex_t IT_0105 = IT_0005*IT_0024;
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0108 = (-2)*IT_0107;
    const ccomplex_t IT_0109 = cpow(s_12 + IT_0012 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0110 = IT_0107*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0108*IT_0111;
    const ccomplex_t IT_0113 = 0.25*IT_0112;
    const ccomplex_t IT_0114 = (-6)*IT_0067;
    const ccomplex_t IT_0115 = IT_0012*IT_0101;
    const ccomplex_t IT_0116 = (-6)*IT_0115;
    const ccomplex_t IT_0117 = 12*IT_0069;
    const ccomplex_t IT_0118 = IT_0114 + IT_0116 + IT_0117;
    const ccomplex_t IT_0119 = m_u*s_14;
    const ccomplex_t IT_0120 = (-18)*IT_0026;
    const ccomplex_t IT_0121 = (-18)*IT_0017;
    const ccomplex_t IT_0122 = (-18)*conj(IT_0026);
    const ccomplex_t IT_0123 = (-18)*conj(IT_0017);
    const ccomplex_t IT_0124 = (-36)*IT_0067;
    const ccomplex_t IT_0125 = (-36)*IT_0115;
    const ccomplex_t IT_0126 = 72*IT_0069;
    const ccomplex_t IT_0127 = IT_0124 + IT_0125 + IT_0126;
    const ccomplex_t IT_0128 = (-0.0833333333333333)*IT_0112;
    const ccomplex_t IT_0129 = cpow(IT_0128, 2);
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0131 = 0.666666666666667*IT_0130;
    const ccomplex_t IT_0132 = (-1.33333333333333)*IT_0130;
    const ccomplex_t IT_0133 = IT_0131*IT_0132;
    const ccomplex_t IT_0134 = IT_0109*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = 0.5*IT_0135;
    const ccomplex_t IT_0137 = (-2)*IT_0128;
    const ccomplex_t IT_0138 = (-12)*IT_0067;
    const ccomplex_t IT_0139 = (-12)*IT_0115;
    const ccomplex_t IT_0140 = 24*IT_0069;
    const ccomplex_t IT_0141 = IT_0138 + IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = -IT_0113;
    const ccomplex_t IT_0143 = (-18)*IT_0115;
    const ccomplex_t IT_0144 = IT_0068 + IT_0070 + IT_0143;
    const ccomplex_t IT_0145 = -IT_0136;
    const ccomplex_t IT_0146 = -s_14;
    const ccomplex_t IT_0147 = s_13 + IT_0146;
    const ccomplex_t IT_0148 = m_u*IT_0147;
    const ccomplex_t IT_0149 = 36*conj(IT_0100);
    const ccomplex_t IT_0150 = 18*IT_0106;
    const ccomplex_t IT_0151 = IT_0100*(18*conj(IT_0017) + 18*conj(IT_0106)) +
       (-6)*(IT_0025 + -conj(IT_0025) + IT_0026 + -conj(IT_0026) + 2*IT_0100 + (
      -2)*conj(IT_0100))*IT_0113 + (IT_0128 + IT_0136)*((-18)*IT_0025 + 18*conj
      (IT_0025) + 18*conj(IT_0026) + (-36)*IT_0100 + IT_0120 + IT_0149) + conj
      (IT_0100)*(18*IT_0017 + IT_0150);
    const ccomplex_t IT_0152 = 18*IT_0101 + IT_0102;
    const ccomplex_t IT_0153 = (-18)*s_34 + IT_0103;
    const ccomplex_t IT_0154 = IT_0100*IT_0104;
    const ccomplex_t IT_0155 = conj(IT_0017)*(18*IT_0026*IT_0027 + (-18)
      *IT_0064*IT_0066 + IT_0017*IT_0071) + conj(IT_0100)*(IT_0100*(36*s_34 + (
      -36)*IT_0101) + (IT_0025 + IT_0026)*IT_0104) + conj(IT_0106)*(18*IT_0025
      *IT_0027 + 18*IT_0064*IT_0066 + IT_0071*IT_0106 + IT_0113*IT_0118 +
       IT_0119*IT_0120 + IT_0115*IT_0121) + IT_0025*(IT_0101*IT_0122 + IT_0119
      *IT_0123) + -IT_0127*(IT_0129 + IT_0136*(IT_0136 + -IT_0137)) + IT_0113*
      (conj(IT_0017)*IT_0118 + -IT_0113*IT_0127 + (-2)*(IT_0136 + (-0.5)*IT_0137
      )*IT_0141) + IT_0017*IT_0118*IT_0142 + IT_0106*(IT_0119*IT_0122 + IT_0115
      *IT_0123 + IT_0118*IT_0142) + -IT_0144*((IT_0017 + -conj(IT_0017) +
       IT_0106 + -conj(IT_0106))*IT_0128 + -(conj(IT_0017) + conj(IT_0106))
      *IT_0136 + -(IT_0017 + IT_0106)*IT_0145) + IT_0148*IT_0151 + conj(IT_0064)
      *(IT_0064*(36*s_34 + 36*IT_0101) + IT_0066*((-18)*IT_0017 + IT_0150) +
       IT_0025*IT_0152 + IT_0026*IT_0153) + conj(IT_0025)*(IT_0025*IT_0102 +
       IT_0101*IT_0120 + IT_0119*IT_0121 + IT_0027*IT_0150 + IT_0064*IT_0152 +
       IT_0154) + conj(IT_0026)*(18*IT_0017*IT_0027 + IT_0026*IT_0102 + IT_0064
      *IT_0153 + IT_0154);
    return create_ccomplex_return(IT_0155);
}

