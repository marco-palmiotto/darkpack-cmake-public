#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_H0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_H0, 2);
    const ccomplex_t IT_0001 = pow(m_sb_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0003 = cpow(V_tb, 2);
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = m_b*mu_h*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = -IT_0009*IT_0012;
    const ccomplex_t IT_0014 = IT_0004*IT_0007;
    const ccomplex_t IT_0015 = IT_0013 + IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -2);
    const ccomplex_t IT_0018 = pow(m_b, 2);
    const ccomplex_t IT_0019 = (-2)*conj(U_sb_00)*(((-0.5)*m_W*IT_0015*
      (IT_0011 + 1./3*IT_0010*IT_0017) + IT_0004*IT_0005*IT_0008*IT_0011*IT_0018
      )*U_sb_01 + -1./2*IT_0006*IT_0008*IT_0009*IT_0011*U_sb_11) + conj(U_sb_10)
      *(IT_0006*IT_0008*IT_0009*IT_0011*U_sb_01 + 0.666666666666667*(m_W*IT_0010
      *IT_0015*IT_0017 + (-3)*IT_0004*IT_0005*IT_0008*IT_0011*IT_0018)*U_sb_11);
    const ccomplex_t IT_0020 = (0 + _Complex_I*-1)*(e_em*IT_0019 + (
      -1.4142135623731)*A_b*IT_0003*IT_0004*(U_sb_01*conj(U_sb_10) + conj
      (U_sb_00)*U_sb_11));
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0002*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0016, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0011
      *IT_0016 + 1./3*IT_0010*IT_0024)*U_sb_00*conj(U_sb_01) + IT_0010*IT_0024
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0026 = (-0.666666666666667)*IT_0025;
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0029 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0015*(IT_0011 + 1.
      /3*IT_0010*IT_0017) + IT_0004*IT_0005*IT_0008*IT_0011*IT_0018)*conj
      (U_sb_00) + -1./2*IT_0006*IT_0008*IT_0009*IT_0011*conj(U_sb_10)) + U_sb_10
      *(IT_0006*IT_0008*IT_0009*IT_0011*conj(U_sb_00) + 0.666666666666667*(m_W
      *IT_0010*IT_0015*IT_0017 + (-3)*IT_0004*IT_0005*IT_0008*IT_0011*IT_0018)
      *conj(U_sb_10));
    const ccomplex_t IT_0030 = (0 + _Complex_I*-1)*(e_em*IT_0029 + (
      -1.4142135623731)*A_b*IT_0003*IT_0004*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0028*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0011
      *IT_0016 + 1./3*IT_0010*IT_0024)*U_sb_00*conj(U_sb_00) + IT_0010*IT_0024
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0035 = (-0.666666666666667)*IT_0034;
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_12 + (-2)*IT_0001 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0038 = m_b*e_em*mu_h*IT_0005*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0039 = IT_0011*IT_0038;
    const ccomplex_t IT_0040 = m_b*e_em*mu_h*IT_0005*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0041 = IT_0011*IT_0040;
    const ccomplex_t IT_0042 = A_b*IT_0003*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0043 = IT_0012*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = A_b*IT_0003*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0046 = IT_0012*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (-0.5)*IT_0039 + 0.5*IT_0041 + (-0.5)*IT_0044 +
       0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0037*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = -e_em*IT_0011*(IT_0007*IT_0009 + -IT_0004
      *IT_0012)*IT_0024;
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = pow(m_Z, 2);
    const ccomplex_t IT_0054 = cpow((-2)*s_12 + (-2)*IT_0001 + IT_0053 + 
      -reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0008*IT_0011*IT_0012;
    const ccomplex_t IT_0056 = IT_0038*IT_0055;
    const ccomplex_t IT_0057 = IT_0040*IT_0055;
    const ccomplex_t IT_0058 = IT_0007*IT_0042;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = IT_0007*IT_0045;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (-0.5)*IT_0056 + 0.5*IT_0057 + 0.5*IT_0059 + (
      -0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0054*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = e_em*IT_0011*(IT_0004*IT_0007 + IT_0009*IT_0012
      )*IT_0024;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0064*IT_0066;
    const ccomplex_t IT_0068 = -IT_0027 + -IT_0036 + IT_0052 + IT_0067;
    const ccomplex_t IT_0069 = pow(m_Z, -2);
    const ccomplex_t IT_0070 = pow(s_34, 2);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0001 + -reg_prop, -1);
    const ccomplex_t IT_0073 = cpow(IT_0016, -4);
    const ccomplex_t IT_0074 = IT_0004*IT_0007*IT_0010*IT_0073;
    const ccomplex_t IT_0075 = IT_0009*IT_0010*IT_0012*IT_0073;
    const ccomplex_t IT_0076 = cpow(IT_0010, -3);
    const ccomplex_t IT_0077 = IT_0004*IT_0007*IT_0076;
    const ccomplex_t IT_0078 = IT_0009*IT_0012*IT_0076;
    const ccomplex_t IT_0079 = IT_0004*IT_0007*IT_0011*IT_0017;
    const ccomplex_t IT_0080 = IT_0009*IT_0011*IT_0012*IT_0017;
    const ccomplex_t IT_0081 = cpow(IT_0010, -2);
    const ccomplex_t IT_0082 = cpow(IT_0017 + IT_0081, -1);
    const ccomplex_t IT_0083 = (IT_0074 + IT_0075 + IT_0077 + IT_0078 + 2
      *IT_0079 + 2*IT_0080)*IT_0082;
    const ccomplex_t IT_0084 = m_W*e_em;
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0072*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (-0.333333333333333)*IT_0034;
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = s_14*s_34*IT_0069;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = s_13 + IT_0093;
    const ccomplex_t IT_0095 = (-3)*IT_0094;
    const ccomplex_t IT_0096 = s_24*s_34*IT_0069;
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = s_23 + IT_0097;
    const ccomplex_t IT_0099 = (-3)*IT_0098;
    const ccomplex_t IT_0100 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0101 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0015*(IT_0011 + 1.
      /3*IT_0010*IT_0017) + IT_0004*IT_0005*IT_0008*IT_0011*IT_0018)*conj
      (U_sb_01) + -1./2*IT_0006*IT_0008*IT_0009*IT_0011*conj(U_sb_11)) + U_sb_10
      *(IT_0006*IT_0008*IT_0009*IT_0011*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0010*IT_0015*IT_0017 + (-3)*IT_0004*IT_0005*IT_0008*IT_0011*IT_0018)
      *conj(U_sb_11));
    const ccomplex_t IT_0102 = (0 + _Complex_I*-1)*(e_em*IT_0101 + (
      -1.4142135623731)*A_b*IT_0003*IT_0004*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0100*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0011
      *IT_0016 + 1./3*IT_0010*IT_0024)*conj(U_sb_00)*U_sb_01 + IT_0010*IT_0024
      *conj(U_sb_10)*U_sb_11);
    const ccomplex_t IT_0107 = (-0.666666666666667)*IT_0106;
    const ccomplex_t IT_0108 = IT_0105*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0110 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0015*(IT_0011 + 1.
      /3*IT_0010*IT_0017) + IT_0004*IT_0005*IT_0008*IT_0011*IT_0018)*conj
      (U_sb_00) + -1./2*IT_0006*IT_0008*IT_0009*IT_0011*conj(U_sb_10)) + U_sb_10
      *(IT_0006*IT_0008*IT_0009*IT_0011*conj(U_sb_00) + 0.666666666666667*(m_W
      *IT_0010*IT_0015*IT_0017 + (-3)*IT_0004*IT_0005*IT_0008*IT_0011*IT_0018)
      *conj(U_sb_10));
    const ccomplex_t IT_0111 = (0 + _Complex_I*-1)*(e_em*IT_0110 + (
      -1.4142135623731)*A_b*IT_0003*IT_0004*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = IT_0109*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0035*IT_0114;
    const ccomplex_t IT_0116 = -IT_0035;
    const ccomplex_t IT_0117 = IT_0089 + IT_0116;
    const ccomplex_t IT_0118 = IT_0088*IT_0117;
    const ccomplex_t IT_0119 = IT_0027 + IT_0036 + IT_0108 + IT_0115 + -IT_0118;
    const ccomplex_t IT_0120 = (-0.333333333333333)*conj(IT_0119);
    const ccomplex_t IT_0121 = pow(s_14, 2);
    const ccomplex_t IT_0122 = IT_0069*IT_0121;
    const ccomplex_t IT_0123 = s_14*s_24*IT_0069;
    const ccomplex_t IT_0124 = (-3)*s_12 + 3*IT_0123;
    const ccomplex_t IT_0125 = (-0.333333333333333)*conj(IT_0068);
    const ccomplex_t IT_0126 = pow(s_24, 2);
    const ccomplex_t IT_0127 = IT_0069*IT_0126;
    const ccomplex_t IT_0128 = (-3)*IT_0068*(conj(IT_0068)*(IT_0000 + -IT_0071
      ) + (-0.333333333333333)*conj(IT_0091)*IT_0095 + IT_0099*IT_0120) + (-3)
      *IT_0091*(conj(IT_0091)*(IT_0001 + -IT_0122) + IT_0120*IT_0124 + IT_0095
      *IT_0125) + (-3)*IT_0119*((-0.333333333333333)*conj(IT_0091)*IT_0124 +
       IT_0099*IT_0125 + conj(IT_0119)*(IT_0001 + -IT_0127));
    return create_ccomplex_return(IT_0128);
}

