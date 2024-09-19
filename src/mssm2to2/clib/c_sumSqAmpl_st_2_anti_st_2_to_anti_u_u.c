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
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = m_t*e_em*mu_h*IT_0002*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = m_t*e_em*mu_h*IT_0002*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0004 + 0.5*IT_0006 + (-0.5)*IT_0010 +
       0.5*IT_0013;
    const ccomplex_t IT_0015 = sin(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = m_u*e_em*IT_0001*IT_0002*IT_0007*IT_0016;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = pow(m_st_2, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + (-2)*IT_0020 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = pow(m_Z, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + (-2)*IT_0020 + IT_0024 + 
      -reg_prop, -1);
    const ccomplex_t IT_0026 = m_u*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = IT_0001*IT_0007*IT_0016;
    const ccomplex_t IT_0029 = IT_0003*IT_0028;
    const ccomplex_t IT_0030 = IT_0005*IT_0028;
    const ccomplex_t IT_0031 = IT_0008*IT_0015;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = IT_0011*IT_0015;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0029 + (-0.5)*IT_0030 + (-0.5)*IT_0032 +
       0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0027*IT_0035;
    const ccomplex_t IT_0037 = IT_0025*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = -IT_0023 + -IT_0038;
    const ccomplex_t IT_0040 = pow(m_u, 2);
    const ccomplex_t IT_0041 = cos(theta_W);
    const ccomplex_t IT_0042 = cpow(IT_0041, -1);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0001
      *IT_0041 + -1./3*IT_0000*IT_0042)*U_st_01*conj(U_st_01) + IT_0000*IT_0042
      *U_st_11*conj(U_st_11));
    const ccomplex_t IT_0044 = 1.33333333333333*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*e_em*IT_0000*IT_0042;
    const ccomplex_t IT_0046 = (-0.666666666666667)*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0020 + -reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0044*IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0000*IT_0042;
    const ccomplex_t IT_0053 = e_em*IT_0052;
    const ccomplex_t IT_0054 = IT_0001*IT_0041;
    const ccomplex_t IT_0055 = e_em*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0053 + (-3)*IT_0055);
    const ccomplex_t IT_0057 = (-0.166666666666667)*IT_0056;
    const ccomplex_t IT_0058 = IT_0047*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0044*IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = s_13 + s_14;
    const ccomplex_t IT_0063 = m_u*IT_0062;
    const ccomplex_t IT_0064 = m_u*IT_0044;
    const ccomplex_t IT_0065 = 0.666666666666667*IT_0043;
    const ccomplex_t IT_0066 = m_u*IT_0065;
    const ccomplex_t IT_0067 = -(IT_0049 + -IT_0059)*(IT_0064 + -IT_0066);
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = 18*s_34;
    const ccomplex_t IT_0070 = 18*IT_0040 + IT_0069;
    const ccomplex_t IT_0071 = (IT_0049 + -IT_0059)*(IT_0064 + -IT_0066);
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = (-18)*IT_0040;
    const ccomplex_t IT_0076 = (-18)*s_34 + IT_0075;
    const ccomplex_t IT_0077 = m_u*s_13;
    const ccomplex_t IT_0078 = s_34*IT_0020;
    const ccomplex_t IT_0079 = (-18)*IT_0078;
    const ccomplex_t IT_0080 = s_13*s_14;
    const ccomplex_t IT_0081 = 36*IT_0080;
    const ccomplex_t IT_0082 = IT_0079 + IT_0081;
    const ccomplex_t IT_0083 = cos(alpha);
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_u*e_em*IT_0001*IT_0002
      *IT_0016*IT_0083;
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = m_t*mu_h*IT_0002;
    const ccomplex_t IT_0087 = sin(alpha);
    const ccomplex_t IT_0088 = IT_0007*IT_0087;
    const ccomplex_t IT_0089 = IT_0015*IT_0083;
    const ccomplex_t IT_0090 = IT_0088 + IT_0089;
    const ccomplex_t IT_0091 = cpow(IT_0041, -2);
    const ccomplex_t IT_0092 = pow(m_t, 2);
    const ccomplex_t IT_0093 = 2*U_st_01*(((-0.5)*m_W*IT_0090*(IT_0001 + -1./3
      *IT_0000*IT_0091) + IT_0001*IT_0002*IT_0016*IT_0083*IT_0092)*conj(U_st_01)
       + 1./2*IT_0001*IT_0016*IT_0086*IT_0087*conj(U_st_11)) + U_st_11*(IT_0001
      *IT_0016*IT_0086*IT_0087*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0000
      *IT_0090*IT_0091 + (-1.5)*IT_0001*IT_0002*IT_0016*IT_0083*IT_0092)*conj
      (U_st_11));
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(e_em*IT_0093 +
       1.4142135623731*A_t*IT_0083*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0085*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0020 + -reg_prop, -1);
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*m_u*e_em*IT_0001*IT_0002
      *IT_0016*IT_0087;
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = IT_0007*IT_0083;
    const ccomplex_t IT_0103 = -IT_0015*IT_0087 + IT_0102;
    const ccomplex_t IT_0104 = (-2)*U_st_01*((IT_0001*IT_0002*IT_0016*IT_0087
      *IT_0092 + 0.5*m_W*(IT_0001 + -1./3*IT_0000*IT_0091)*IT_0103)*conj(U_st_01
      ) + -1./2*IT_0001*IT_0016*IT_0083*IT_0086*conj(U_st_11)) + U_st_11*
      (IT_0001*IT_0016*IT_0083*IT_0086*conj(U_st_01) + (-1.33333333333333)*(1.5
      *IT_0001*IT_0002*IT_0016*IT_0087*IT_0092 + m_W*IT_0000*IT_0091*IT_0103)
      *conj(U_st_11));
    const ccomplex_t IT_0105 = (0 + _Complex_I*-1)*(e_em*IT_0104 + (
      -1.4142135623731)*A_t*IT_0087*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = IT_0101*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_12 + (-2)*IT_0020 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = -IT_0099 + -IT_0110;
    const ccomplex_t IT_0112 = IT_0069 + IT_0075;
    const ccomplex_t IT_0113 = -s_14;
    const ccomplex_t IT_0114 = s_13 + IT_0113;
    const ccomplex_t IT_0115 = m_u*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0117 = (-2)*IT_0116;
    const ccomplex_t IT_0118 = cpow(s_12 + IT_0020 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0119 = IT_0116*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = IT_0117*IT_0120;
    const ccomplex_t IT_0122 = 0.25*IT_0121;
    const ccomplex_t IT_0123 = (-0.0833333333333333)*IT_0121;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0125 = 0.666666666666667*IT_0124;
    const ccomplex_t IT_0126 = (-1.33333333333333)*IT_0124;
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = IT_0118*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = 36*conj(IT_0111);
    const ccomplex_t IT_0132 = (-18)*IT_0074;
    const ccomplex_t IT_0133 = (18*conj(IT_0051) + 18*conj(IT_0061))*IT_0111 +
       (18*IT_0051 + 18*IT_0061)*conj(IT_0111) + (-6)*(IT_0068 + -conj(IT_0068) 
      + IT_0074 + -conj(IT_0074) + 2*IT_0111 + (-2)*conj(IT_0111))*IT_0122 + 
      (IT_0123 + IT_0130)*((-18)*IT_0068 + 18*conj(IT_0068) + 18*conj(IT_0074) +
       (-36)*IT_0111 + IT_0131 + IT_0132);
    const ccomplex_t IT_0134 = (-6)*IT_0078;
    const ccomplex_t IT_0135 = IT_0020*IT_0040;
    const ccomplex_t IT_0136 = (-6)*IT_0135;
    const ccomplex_t IT_0137 = 12*IT_0080;
    const ccomplex_t IT_0138 = IT_0134 + IT_0136 + IT_0137;
    const ccomplex_t IT_0139 = m_u*s_14;
    const ccomplex_t IT_0140 = (-18)*IT_0061;
    const ccomplex_t IT_0141 = (-18)*conj(IT_0074);
    const ccomplex_t IT_0142 = (-18)*conj(IT_0061);
    const ccomplex_t IT_0143 = (-36)*IT_0078;
    const ccomplex_t IT_0144 = (-36)*IT_0135;
    const ccomplex_t IT_0145 = 72*IT_0080;
    const ccomplex_t IT_0146 = IT_0143 + IT_0144 + IT_0145;
    const ccomplex_t IT_0147 = cpow(IT_0123, 2);
    const ccomplex_t IT_0148 = (-2)*IT_0123;
    const ccomplex_t IT_0149 = (-12)*IT_0078;
    const ccomplex_t IT_0150 = (-12)*IT_0135;
    const ccomplex_t IT_0151 = 24*IT_0080;
    const ccomplex_t IT_0152 = IT_0149 + IT_0150 + IT_0151;
    const ccomplex_t IT_0153 = -IT_0122;
    const ccomplex_t IT_0154 = (-18)*IT_0135;
    const ccomplex_t IT_0155 = IT_0079 + IT_0081 + IT_0154;
    const ccomplex_t IT_0156 = -IT_0130;
    const ccomplex_t IT_0157 = IT_0111*IT_0112;
    const ccomplex_t IT_0158 = conj(IT_0039)*(IT_0039*(36*s_34 + 36*IT_0040) +
       (18*IT_0051 + (-18)*IT_0061)*IT_0063 + IT_0068*IT_0070 + IT_0074*IT_0076)
       + conj(IT_0061)*((-18)*IT_0039*IT_0063 + 18*IT_0074*IT_0077 + IT_0061
      *IT_0082) + conj(IT_0111)*((36*s_34 + (-36)*IT_0040)*IT_0111 + (IT_0068 +
       IT_0074)*IT_0112) + IT_0115*IT_0133 + conj(IT_0051)*(18*IT_0039*IT_0063 +
       18*IT_0068*IT_0077 + IT_0051*IT_0082 + IT_0122*IT_0138 + IT_0132*IT_0139 
      + IT_0135*IT_0140) + IT_0068*(IT_0040*IT_0141 + IT_0139*IT_0142) + 
      -IT_0146*(IT_0147 + IT_0130*(IT_0130 + -IT_0148)) + IT_0122*(conj(IT_0061)
      *IT_0138 + -IT_0122*IT_0146 + (-2)*(IT_0130 + (-0.5)*IT_0148)*IT_0152) +
       IT_0061*IT_0138*IT_0153 + IT_0051*(IT_0139*IT_0141 + IT_0135*IT_0142 +
       IT_0138*IT_0153) + -IT_0155*((IT_0051 + -conj(IT_0051) + IT_0061 + -conj
      (IT_0061))*IT_0123 + -(conj(IT_0051) + conj(IT_0061))*IT_0130 + -(IT_0051 
      + IT_0061)*IT_0156) + conj(IT_0074)*(IT_0069*IT_0074 + IT_0039*IT_0076 +
       18*IT_0061*IT_0077 + IT_0157) + conj(IT_0068)*(IT_0068*IT_0069 + IT_0039
      *IT_0070 + 18*IT_0051*IT_0077 + IT_0040*IT_0132 + IT_0139*IT_0140 +
       IT_0157);
    return create_ccomplex_return(IT_0158);
}

