#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_A0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
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
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_sb_1, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0002 = cpow(V_tb, 2);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_b*mu_h*IT_0004;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = sin(beta);
    const ccomplex_t IT_0012 = -IT_0008*IT_0011;
    const ccomplex_t IT_0013 = IT_0003*IT_0006;
    const ccomplex_t IT_0014 = IT_0012 + IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -2);
    const ccomplex_t IT_0017 = pow(m_b, 2);
    const ccomplex_t IT_0018 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0014*(IT_0010 + 1.
      /3*IT_0009*IT_0016) + IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)*conj
      (U_sb_00) + -1./2*IT_0005*IT_0007*IT_0008*IT_0010*conj(U_sb_10)) + U_sb_10
      *(IT_0005*IT_0007*IT_0008*IT_0010*conj(U_sb_00) + 0.666666666666667*(m_W
      *IT_0009*IT_0014*IT_0016 + (-3)*IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)
      *conj(U_sb_10));
    const ccomplex_t IT_0019 = (0 + _Complex_I*-1)*(e_em*IT_0018 + (
      -1.4142135623731)*A_b*IT_0002*IT_0003*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0001*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = cpow(IT_0015, -1);
    const ccomplex_t IT_0024 = -e_em*IT_0010*(IT_0006*IT_0008 + -IT_0003
      *IT_0011)*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = pow(m_A0, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0026 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0028 = m_b*e_em*mu_h*IT_0004*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0029 = IT_0010*IT_0028;
    const ccomplex_t IT_0030 = m_b*e_em*mu_h*IT_0004*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0031 = IT_0010*IT_0030;
    const ccomplex_t IT_0032 = A_b*IT_0002*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0033 = IT_0011*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = A_b*IT_0002*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0036 = IT_0011*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (-0.5)*IT_0029 + 0.5*IT_0031 + (-0.5)*IT_0034 +
       0.5*IT_0037;
    const ccomplex_t IT_0039 = IT_0027*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0010
      *IT_0015 + 1./3*IT_0009*IT_0023)*U_sb_00*conj(U_sb_00) + IT_0009*IT_0023
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0042 = (-0.666666666666667)*IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0006*IT_0008;
    const ccomplex_t IT_0046 = IT_0003*IT_0011;
    const ccomplex_t IT_0047 = IT_0045 + IT_0046;
    const ccomplex_t IT_0048 = 2*U_sb_00*((IT_0004*IT_0007*IT_0008*IT_0010
      *IT_0017 + (-0.5)*m_W*(IT_0010 + 1./3*IT_0009*IT_0016)*IT_0047)*conj
      (U_sb_00) + 1./2*IT_0003*IT_0005*IT_0007*IT_0010*conj(U_sb_10)) + U_sb_10*
      (IT_0003*IT_0005*IT_0007*IT_0010*conj(U_sb_00) + 2*(IT_0004*IT_0007
      *IT_0008*IT_0010*IT_0017 + -1./3*m_W*IT_0009*IT_0016*IT_0047)*conj(U_sb_10
      ));
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(e_em*IT_0048 +
       1.4142135623731*A_b*IT_0002*IT_0008*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = IT_0044*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = e_em*IT_0010*(IT_0003*IT_0006 + IT_0008*IT_0011
      )*IT_0023;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0000 + IT_0026 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0057 = m_b*e_em*mu_h*IT_0004*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0058 = IT_0010*IT_0057;
    const ccomplex_t IT_0059 = m_b*e_em*mu_h*IT_0004*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0060 = IT_0010*IT_0059;
    const ccomplex_t IT_0061 = A_b*IT_0002*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0062 = IT_0011*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = A_b*IT_0002*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0065 = IT_0011*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0058 + (-0.5)*IT_0060 + 0.5*IT_0063 + (
      -0.5)*IT_0066;
    const ccomplex_t IT_0068 = IT_0056*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0010
      *IT_0015 + 1./3*IT_0009*IT_0023)*U_sb_00*conj(U_sb_01) + IT_0009*IT_0023
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0071 = (-0.666666666666667)*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_23 + IT_0026 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0038*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0042*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_13 + IT_0000 + IT_0026 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0078 = m_b*e_em*mu_h*IT_0004*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0079 = IT_0010*IT_0078;
    const ccomplex_t IT_0080 = m_b*e_em*mu_h*IT_0004*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0081 = IT_0010*IT_0080;
    const ccomplex_t IT_0082 = A_b*IT_0002*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0083 = IT_0011*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = A_b*IT_0002*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0086 = IT_0011*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*IT_0079 + 0.5*IT_0081 + (-0.5)*IT_0084 +
       0.5*IT_0087;
    const ccomplex_t IT_0089 = IT_0077*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0010
      *IT_0015 + 1./3*IT_0009*IT_0023)*conj(U_sb_00)*U_sb_01 + IT_0009*IT_0023
      *conj(U_sb_10)*U_sb_11);
    const ccomplex_t IT_0092 = (-0.666666666666667)*IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = -IT_0025 + IT_0043 + -IT_0055 + IT_0072 +
       IT_0076 + IT_0093;
    const ccomplex_t IT_0095 = pow(m_Z, -2);
    const ccomplex_t IT_0096 = pow(s_24, 2);
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = -IT_0025 + -IT_0055;
    const ccomplex_t IT_0099 = s_14*s_24*IT_0095;
    const ccomplex_t IT_0100 = (-3)*s_12 + 3*IT_0099;
    const ccomplex_t IT_0101 = s_24*s_34*IT_0095;
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = s_23 + IT_0102;
    const ccomplex_t IT_0104 = (-3)*IT_0103;
    const ccomplex_t IT_0105 = -IT_0072 + -IT_0076;
    const ccomplex_t IT_0106 = (-0.333333333333333)*conj(IT_0105);
    const ccomplex_t IT_0107 = pow(s_14, 2);
    const ccomplex_t IT_0108 = IT_0095*IT_0107;
    const ccomplex_t IT_0109 = s_14*s_34*IT_0095;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = s_13 + IT_0110;
    const ccomplex_t IT_0112 = (-3)*IT_0111;
    const ccomplex_t IT_0113 = (-0.333333333333333)*conj(IT_0094);
    const ccomplex_t IT_0114 = pow(s_34, 2);
    const ccomplex_t IT_0115 = IT_0095*IT_0114;
    const ccomplex_t IT_0116 = (-3)*IT_0094*(conj(IT_0094)*(IT_0000 + -IT_0097
      ) + (-0.333333333333333)*conj(IT_0098)*IT_0100 + IT_0104*IT_0106) + (-3)
      *IT_0098*(conj(IT_0098)*(IT_0000 + -IT_0108) + IT_0106*IT_0112 + IT_0100
      *IT_0113) + (-3)*IT_0105*((-0.333333333333333)*conj(IT_0098)*IT_0112 +
       IT_0104*IT_0113 + conj(IT_0105)*(IT_0026 + -IT_0115));
    return create_ccomplex_return(IT_0116);
}

