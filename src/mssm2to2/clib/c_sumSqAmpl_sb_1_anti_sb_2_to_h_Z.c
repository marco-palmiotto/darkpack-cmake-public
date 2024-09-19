#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_h_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_h_Z(
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_00*conj(U_sb_01) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0005 = (-0.333333333333333)*IT_0004;
    const ccomplex_t IT_0006 = m_W*e_em;
    const ccomplex_t IT_0007 = cos(alpha);
    const ccomplex_t IT_0008 = cpow(IT_0000, -4);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = IT_0002*IT_0007*IT_0008*IT_0009;
    const ccomplex_t IT_0011 = cos(beta);
    const ccomplex_t IT_0012 = sin(alpha);
    const ccomplex_t IT_0013 = IT_0002*IT_0008*IT_0011*IT_0012;
    const ccomplex_t IT_0014 = cpow(IT_0002, -3);
    const ccomplex_t IT_0015 = IT_0007*IT_0009*IT_0014;
    const ccomplex_t IT_0016 = IT_0011*IT_0012*IT_0014;
    const ccomplex_t IT_0017 = cpow(IT_0000, -2);
    const ccomplex_t IT_0018 = IT_0003*IT_0007*IT_0009*IT_0017;
    const ccomplex_t IT_0019 = IT_0003*IT_0011*IT_0012*IT_0017;
    const ccomplex_t IT_0020 = cpow(IT_0002, -2);
    const ccomplex_t IT_0021 = cpow(IT_0017 + IT_0020, -1);
    const ccomplex_t IT_0022 = (IT_0010 + -IT_0013 + IT_0015 + -IT_0016 + 2
      *IT_0018 + (-2)*IT_0019)*IT_0021;
    const ccomplex_t IT_0023 = IT_0006*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = pow(m_sb_1, 2);
    const ccomplex_t IT_0026 = pow(m_sb_2, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0025 + -IT_0026 + -reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0024*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0005*IT_0029;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = pow(m_Z, -2);
    const ccomplex_t IT_0033 = pow(s_14, 2);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_00*conj(U_sb_00) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0036 = (-0.666666666666667)*IT_0035;
    const ccomplex_t IT_0037 = cpow(V_tb, 2);
    const ccomplex_t IT_0038 = pow(m_W, -1);
    const ccomplex_t IT_0039 = m_b*mu_h*IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0011, -1);
    const ccomplex_t IT_0041 = IT_0011*IT_0012;
    const ccomplex_t IT_0042 = IT_0007*IT_0009;
    const ccomplex_t IT_0043 = IT_0041 + IT_0042;
    const ccomplex_t IT_0044 = pow(m_b, 2);
    const ccomplex_t IT_0045 = 2*U_sb_00*(((-0.5)*m_W*(IT_0003 + 1./3*IT_0002
      *IT_0017)*IT_0043 + IT_0003*IT_0012*IT_0038*IT_0040*IT_0044)*conj(U_sb_01)
       + 1./2*IT_0003*IT_0007*IT_0039*IT_0040*conj(U_sb_11)) + U_sb_10*(IT_0003
      *IT_0007*IT_0039*IT_0040*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0002
      *IT_0017*IT_0043 + (-3)*IT_0003*IT_0012*IT_0038*IT_0040*IT_0044)*conj
      (U_sb_11));
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(e_em*IT_0045 +
       1.4142135623731*A_b*IT_0012*IT_0037*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = pow(m_h, 2);
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0026 + IT_0048 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0036*IT_0051;
    const ccomplex_t IT_0053 = e_em*IT_0001*IT_0003*(IT_0007*IT_0011 + IT_0009
      *IT_0012);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = m_b*e_em*mu_h*IT_0038*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0056 = IT_0003*IT_0055;
    const ccomplex_t IT_0057 = m_b*e_em*mu_h*IT_0038*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0058 = IT_0003*IT_0057;
    const ccomplex_t IT_0059 = A_b*IT_0037*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0060 = IT_0009*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = A_b*IT_0037*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0063 = IT_0009*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (-0.5)*IT_0056 + 0.5*IT_0058 + (-0.5)*IT_0061 +
       0.5*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_12 + -IT_0025 + -IT_0026 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0054*IT_0068;
    const ccomplex_t IT_0070 = e_em*IT_0001*IT_0003*(IT_0007*IT_0009 + 
      -IT_0011*IT_0012);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0003*IT_0009*IT_0040;
    const ccomplex_t IT_0073 = IT_0055*IT_0072;
    const ccomplex_t IT_0074 = IT_0057*IT_0072;
    const ccomplex_t IT_0075 = IT_0011*IT_0059;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = IT_0011*IT_0062;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0073 + 0.5*IT_0074 + 0.5*IT_0076 + (
      -0.5)*IT_0078;
    const ccomplex_t IT_0080 = pow(m_Z, 2);
    const ccomplex_t IT_0081 = cpow((-2)*s_12 + -IT_0025 + -IT_0026 + IT_0080 
      + -reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = IT_0071*IT_0083;
    const ccomplex_t IT_0085 = (-0.666666666666667)*IT_0004;
    const ccomplex_t IT_0086 = 2*U_sb_01*(((-0.5)*m_W*(IT_0003 + 1./3*IT_0002
      *IT_0017)*IT_0043 + IT_0003*IT_0012*IT_0038*IT_0040*IT_0044)*conj(U_sb_01)
       + 1./2*IT_0003*IT_0007*IT_0039*IT_0040*conj(U_sb_11)) + U_sb_11*(IT_0003
      *IT_0007*IT_0039*IT_0040*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0002
      *IT_0017*IT_0043 + (-3)*IT_0003*IT_0012*IT_0038*IT_0040*IT_0044)*conj
      (U_sb_11));
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(e_em*IT_0086 +
       1.4142135623731*A_b*IT_0012*IT_0037*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_23 + IT_0048 + (0 + _Complex_I*1)
      *m_sb_2*Gamma_b2 + reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0085*IT_0091;
    const ccomplex_t IT_0093 = -IT_0052 + IT_0069 + IT_0084 + -IT_0092;
    const ccomplex_t IT_0094 = s_14*s_34*IT_0032;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = s_13 + IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_01*conj(U_sb_01) + IT_0001*IT_0002
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0098 = (-0.666666666666667)*IT_0097;
    const ccomplex_t IT_0099 = 2*U_sb_00*(((-0.5)*m_W*(IT_0003 + 1./3*IT_0002
      *IT_0017)*IT_0043 + IT_0003*IT_0012*IT_0038*IT_0040*IT_0044)*conj(U_sb_01)
       + 1./2*IT_0003*IT_0007*IT_0039*IT_0040*conj(U_sb_11)) + U_sb_10*(IT_0003
      *IT_0007*IT_0039*IT_0040*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0002
      *IT_0017*IT_0043 + (-3)*IT_0003*IT_0012*IT_0038*IT_0040*IT_0044)*conj
      (U_sb_11));
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(e_em*IT_0099 +
       1.4142135623731*A_b*IT_0012*IT_0037*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_13 + IT_0025 + IT_0048 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0098*IT_0104;
    const ccomplex_t IT_0106 = 2*U_sb_00*(((-0.5)*m_W*(IT_0003 + 1./3*IT_0002
      *IT_0017)*IT_0043 + IT_0003*IT_0012*IT_0038*IT_0040*IT_0044)*conj(U_sb_00)
       + 1./2*IT_0003*IT_0007*IT_0039*IT_0040*conj(U_sb_10)) + U_sb_10*(IT_0003
      *IT_0007*IT_0039*IT_0040*conj(U_sb_00) + (-0.666666666666667)*(m_W*IT_0002
      *IT_0017*IT_0043 + (-3)*IT_0003*IT_0012*IT_0038*IT_0040*IT_0044)*conj
      (U_sb_10));
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(e_em*IT_0106 +
       1.4142135623731*A_b*IT_0012*IT_0037*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_13 + IT_0048 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0110 = IT_0108*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0085*IT_0111;
    const ccomplex_t IT_0113 = -IT_0085;
    const ccomplex_t IT_0114 = IT_0005 + IT_0113;
    const ccomplex_t IT_0115 = IT_0029*IT_0114;
    const ccomplex_t IT_0116 = IT_0052 + IT_0092 + IT_0105 + IT_0112 + -IT_0115;
    const ccomplex_t IT_0117 = s_14*s_24*IT_0032;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = s_12 + IT_0118;
    const ccomplex_t IT_0120 = (-3)*conj(IT_0031);
    const ccomplex_t IT_0121 = pow(s_34, 2);
    const ccomplex_t IT_0122 = IT_0032*IT_0121;
    const ccomplex_t IT_0123 = s_24*s_34*IT_0032;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = s_23 + IT_0124;
    const ccomplex_t IT_0126 = pow(s_24, 2);
    const ccomplex_t IT_0127 = IT_0032*IT_0126;
    const ccomplex_t IT_0128 = (-3)*IT_0116;
    const ccomplex_t IT_0129 = (-3)*IT_0031*(conj(IT_0031)*(IT_0025 + -IT_0034
      ) + conj(IT_0093)*IT_0096 + conj(IT_0116)*IT_0119) + IT_0116*IT_0119
      *IT_0120 + (-3)*IT_0093*((-0.333333333333333)*IT_0096*IT_0120 + conj
      (IT_0093)*(IT_0048 + -IT_0122) + conj(IT_0116)*IT_0125) + (conj(IT_0093)
      *IT_0125 + conj(IT_0116)*(IT_0026 + -IT_0127))*IT_0128;
    return create_ccomplex_return(IT_0129);
}

