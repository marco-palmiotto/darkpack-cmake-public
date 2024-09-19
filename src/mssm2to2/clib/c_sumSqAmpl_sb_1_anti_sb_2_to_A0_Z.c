#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_A0_Z(
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
    const ccomplex_t IT_0000 = pow(m_A0, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_b*e_em*mu_h*IT_0004*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = m_b*e_em*mu_h*IT_0004*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(V_tb, 2);
    const ccomplex_t IT_0011 = A_b*IT_0010*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = A_b*IT_0010*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0006 + 0.5*IT_0008 + (-0.5)*IT_0013 +
       0.5*IT_0016;
    const ccomplex_t IT_0018 = IT_0001*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cos(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0003
      *IT_0020 + 1./3*IT_0002*IT_0021)*U_sb_00*conj(U_sb_01) + IT_0002*IT_0021
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0023 = (-0.666666666666667)*IT_0022;
    const ccomplex_t IT_0024 = IT_0019*IT_0023;
    const ccomplex_t IT_0025 = pow(m_sb_2, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0000 + IT_0025 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0027 = m_b*e_em*mu_h*IT_0004*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = m_b*e_em*mu_h*IT_0004*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0030 = IT_0003*IT_0029;
    const ccomplex_t IT_0031 = A_b*IT_0010*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0032 = IT_0009*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = A_b*IT_0010*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0035 = IT_0009*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0028 + 0.5*IT_0030 + (-0.5)*IT_0033 +
       0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0026*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0003
      *IT_0020 + 1./3*IT_0002*IT_0021)*U_sb_00*conj(U_sb_00) + IT_0002*IT_0021
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0041 = (-0.666666666666667)*IT_0040;
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = -IT_0024 + -IT_0042;
    const ccomplex_t IT_0044 = pow(m_Z, -2);
    const ccomplex_t IT_0045 = pow(s_34, 2);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = pow(m_sb_1, 2);
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + -IT_0025 + -IT_0047 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0049 = cos(alpha);
    const ccomplex_t IT_0050 = m_b*mu_h*IT_0004;
    const ccomplex_t IT_0051 = cos(beta);
    const ccomplex_t IT_0052 = cpow(IT_0051, -1);
    const ccomplex_t IT_0053 = sin(alpha);
    const ccomplex_t IT_0054 = -IT_0009*IT_0053;
    const ccomplex_t IT_0055 = IT_0049*IT_0051;
    const ccomplex_t IT_0056 = IT_0054 + IT_0055;
    const ccomplex_t IT_0057 = cpow(IT_0020, -2);
    const ccomplex_t IT_0058 = pow(m_b, 2);
    const ccomplex_t IT_0059 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0056*(IT_0003 + 1.
      /3*IT_0002*IT_0057) + IT_0003*IT_0004*IT_0049*IT_0052*IT_0058)*conj
      (U_sb_01) + -1./2*IT_0003*IT_0050*IT_0052*IT_0053*conj(U_sb_11)) + U_sb_10
      *(IT_0003*IT_0050*IT_0052*IT_0053*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0002*IT_0056*IT_0057 + (-3)*IT_0003*IT_0004*IT_0049*IT_0052*IT_0058)
      *conj(U_sb_11));
    const ccomplex_t IT_0060 = (0 + _Complex_I*-1)*(e_em*IT_0059 + (
      -1.4142135623731)*A_b*IT_0010*IT_0049*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0048*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = e_em*IT_0003*IT_0021*(IT_0009*IT_0049 + 
      -IT_0051*IT_0053);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0025 + -IT_0047 + -reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0051*IT_0053;
    const ccomplex_t IT_0068 = IT_0009*IT_0049;
    const ccomplex_t IT_0069 = IT_0067 + IT_0068;
    const ccomplex_t IT_0070 = 2*U_sb_00*((IT_0003*IT_0004*IT_0052*IT_0053
      *IT_0058 + (-0.5)*m_W*(IT_0003 + 1./3*IT_0002*IT_0057)*IT_0069)*conj
      (U_sb_01) + 1./2*IT_0003*IT_0049*IT_0050*IT_0052*conj(U_sb_11)) + U_sb_10*
      (IT_0003*IT_0049*IT_0050*IT_0052*conj(U_sb_01) + 2*(IT_0003*IT_0004
      *IT_0052*IT_0053*IT_0058 + -1./3*m_W*IT_0002*IT_0057*IT_0069)*conj(U_sb_11
      ));
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(e_em*IT_0070 +
       1.4142135623731*A_b*IT_0010*IT_0053*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = IT_0066*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = e_em*IT_0003*IT_0021*(IT_0049*IT_0051 + IT_0009
      *IT_0053);
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = -IT_0065 + -IT_0077;
    const ccomplex_t IT_0079 = s_14*s_34*IT_0044;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = s_13 + IT_0080;
    const ccomplex_t IT_0082 = (-3)*IT_0081;
    const ccomplex_t IT_0083 = s_24*s_34*IT_0044;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = s_23 + IT_0084;
    const ccomplex_t IT_0086 = (-3)*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0088 = m_b*e_em*mu_h*IT_0004*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0089 = IT_0003*IT_0088;
    const ccomplex_t IT_0090 = m_b*e_em*mu_h*IT_0004*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0091 = IT_0003*IT_0090;
    const ccomplex_t IT_0092 = A_b*IT_0010*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0093 = IT_0009*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = A_b*IT_0010*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0096 = IT_0009*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0089 + 0.5*IT_0091 + (-0.5)*IT_0094 +
       0.5*IT_0097;
    const ccomplex_t IT_0099 = IT_0087*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0023*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_13 + IT_0000 + IT_0047 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0037*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0003
      *IT_0020 + 1./3*IT_0002*IT_0021)*U_sb_01*conj(U_sb_01) + IT_0002*IT_0021
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0106 = (-0.666666666666667)*IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = IT_0024 + IT_0042 + -IT_0065 + -IT_0077 +
       IT_0101 + IT_0107;
    const ccomplex_t IT_0109 = (-0.333333333333333)*conj(IT_0108);
    const ccomplex_t IT_0110 = pow(s_14, 2);
    const ccomplex_t IT_0111 = IT_0044*IT_0110;
    const ccomplex_t IT_0112 = s_14*s_24*IT_0044;
    const ccomplex_t IT_0113 = (-3)*s_12 + 3*IT_0112;
    const ccomplex_t IT_0114 = (-0.333333333333333)*conj(IT_0043);
    const ccomplex_t IT_0115 = pow(s_24, 2);
    const ccomplex_t IT_0116 = IT_0044*IT_0115;
    const ccomplex_t IT_0117 = (-3)*IT_0043*(conj(IT_0043)*(IT_0000 + -IT_0046
      ) + (-0.333333333333333)*conj(IT_0078)*IT_0082 + IT_0086*IT_0109) + (-3)
      *IT_0078*(conj(IT_0078)*(IT_0047 + -IT_0111) + IT_0109*IT_0113 + IT_0082
      *IT_0114) + (-3)*IT_0108*((-0.333333333333333)*conj(IT_0078)*IT_0113 +
       IT_0086*IT_0114 + conj(IT_0108)*(IT_0025 + -IT_0116));
    return create_ccomplex_return(IT_0117);
}

