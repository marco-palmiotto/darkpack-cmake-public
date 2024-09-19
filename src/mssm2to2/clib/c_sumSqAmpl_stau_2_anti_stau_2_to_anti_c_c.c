#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_stau_2_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_stau_2_to_anti_c_c(
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
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*U_stau_01*conj(U_stau_01) + IT_0001*IT_0002*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0005 = (-2)*IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0007 = (-0.666666666666667)*IT_0006;
    const ccomplex_t IT_0008 = pow(m_stau_2, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = pow(m_c, 2);
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = m_c*s_14;
    const ccomplex_t IT_0017 = -IT_0004;
    const ccomplex_t IT_0018 = m_c*IT_0017;
    const ccomplex_t IT_0019 = m_c*IT_0005;
    const ccomplex_t IT_0020 = IT_0001*IT_0002;
    const ccomplex_t IT_0021 = e_em*IT_0020;
    const ccomplex_t IT_0022 = IT_0000*IT_0003;
    const ccomplex_t IT_0023 = e_em*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0021 + (-3)*IT_0023);
    const ccomplex_t IT_0025 = (-0.166666666666667)*IT_0024;
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -(IT_0018 + -IT_0019)*(IT_0011 + -IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = IT_0005*IT_0027;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = (-6)*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_12 + (-2)*IT_0008 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0035 = pow(m_W, -1);
    const ccomplex_t IT_0036 = sin(beta);
    const ccomplex_t IT_0037 = cpow(IT_0036, -1);
    const ccomplex_t IT_0038 = sin(alpha);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*m_c*e_em*IT_0003*IT_0035
      *IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = cos(alpha);
    const ccomplex_t IT_0042 = mu_h*m_tau*IT_0035;
    const ccomplex_t IT_0043 = cos(beta);
    const ccomplex_t IT_0044 = cpow(IT_0043, -1);
    const ccomplex_t IT_0045 = -IT_0036*IT_0038;
    const ccomplex_t IT_0046 = IT_0041*IT_0043;
    const ccomplex_t IT_0047 = IT_0045 + IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0000, -2);
    const ccomplex_t IT_0049 = pow(m_tau, 2);
    const ccomplex_t IT_0050 = (-2)*U_stau_01*(((-0.5)*m_W*IT_0047*(IT_0003 + 
      -IT_0002*IT_0048) + IT_0003*IT_0035*IT_0041*IT_0044*IT_0049)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0038*IT_0042*IT_0044*conj(U_stau_11)) +
       U_stau_11*(IT_0003*IT_0038*IT_0042*IT_0044*conj(U_stau_01) + 2*(m_W
      *IT_0002*IT_0047*IT_0048 + -IT_0003*IT_0035*IT_0041*IT_0044*IT_0049)*conj
      (U_stau_11));
    const ccomplex_t IT_0051 = (0 + _Complex_I*-1)*(e_em*IT_0050 + (
      -1.4142135623731)*A_tau*IT_0041*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = IT_0040*IT_0052;
    const ccomplex_t IT_0054 = IT_0034*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_c*e_em*IT_0003*IT_0035
      *IT_0037*IT_0041;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0038*IT_0043;
    const ccomplex_t IT_0060 = IT_0036*IT_0041;
    const ccomplex_t IT_0061 = IT_0059 + IT_0060;
    const ccomplex_t IT_0062 = 2*U_stau_01*((IT_0003*IT_0035*IT_0038*IT_0044
      *IT_0049 + (-0.5)*m_W*(IT_0003 + -IT_0002*IT_0048)*IT_0061)*conj(U_stau_01
      ) + 1./2*IT_0003*IT_0041*IT_0042*IT_0044*conj(U_stau_11)) + U_stau_11*
      (IT_0003*IT_0041*IT_0042*IT_0044*conj(U_stau_01) + 2*(IT_0003*IT_0035
      *IT_0038*IT_0044*IT_0049 + -m_W*IT_0002*IT_0048*IT_0061)*conj(U_stau_11));
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(e_em*IT_0062 +
       1.4142135623731*A_tau*IT_0038*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0058*IT_0064;
    const ccomplex_t IT_0066 = IT_0056*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0055 + -IT_0067;
    const ccomplex_t IT_0069 = -s_14;
    const ccomplex_t IT_0070 = s_13 + IT_0069;
    const ccomplex_t IT_0071 = m_c*IT_0070;
    const ccomplex_t IT_0072 = -IT_0030;
    const ccomplex_t IT_0073 = -IT_0014;
    const ccomplex_t IT_0074 = s_34 + IT_0073;
    const ccomplex_t IT_0075 = 12*s_34;
    const ccomplex_t IT_0076 = m_c*e_em*IT_0003*IT_0035;
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = e_em*mu_h*m_tau*IT_0035*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0079 = IT_0003*IT_0036*IT_0044;
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = e_em*mu_h*m_tau*IT_0035*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0084 = IT_0043*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0087 = IT_0043*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (-0.5)*IT_0080 + 0.5*IT_0082 + 0.5*IT_0085 + (
      -0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0077*IT_0089;
    const ccomplex_t IT_0091 = pow(m_Z, 2);
    const ccomplex_t IT_0092 = cpow((-2)*s_12 + (-2)*IT_0008 + IT_0091 + 
      -reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = IT_0003*IT_0078;
    const ccomplex_t IT_0096 = IT_0003*IT_0081;
    const ccomplex_t IT_0097 = IT_0036*IT_0083;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = IT_0036*IT_0086;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = (-0.5)*IT_0095 + 0.5*IT_0096 + (-0.5)*IT_0098 +
       0.5*IT_0100;
    const ccomplex_t IT_0102 = m_c*e_em*IT_0003*IT_0035*IT_0037*IT_0043;
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0101*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_12 + (-2)*IT_0008 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = -IT_0094 + -IT_0107;
    const ccomplex_t IT_0109 = s_34 + IT_0014;
    const ccomplex_t IT_0110 = (-6)*conj(IT_0072);
    const ccomplex_t IT_0111 = s_13 + s_14;
    const ccomplex_t IT_0112 = m_c*IT_0111;
    const ccomplex_t IT_0113 = -conj(IT_0032);
    const ccomplex_t IT_0114 = cpow(s_12 + IT_0008 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0116 = 2*IT_0115;
    const ccomplex_t IT_0117 = 0.666666666666667*IT_0115;
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = IT_0114*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = s_13*s_14;
    const ccomplex_t IT_0123 = s_34*IT_0008;
    const ccomplex_t IT_0124 = (-6)*IT_0072;
    const ccomplex_t IT_0125 = (-6)*IT_0123;
    const ccomplex_t IT_0126 = (-6)*IT_0015;
    const ccomplex_t IT_0127 = 12*IT_0122;
    const ccomplex_t IT_0128 = IT_0125 + IT_0126 + IT_0127;
    const ccomplex_t IT_0129 = IT_0121*((-12)*IT_0015 + 24*IT_0122 + (-12)
      *IT_0123) + 6*IT_0071*(IT_0030 + -conj(IT_0030) + 2*IT_0068 + (-2)*conj
      (IT_0068) + -conj(IT_0072) + (-0.166666666666667)*IT_0124) + (IT_0013 + 
      -conj(IT_0013) + IT_0032 + IT_0113)*IT_0128;
    const ccomplex_t IT_0130 = conj(IT_0108)*IT_0109;
    const ccomplex_t IT_0131 = m_c*s_13;
    const ccomplex_t IT_0132 = IT_0125 + IT_0127;
    const ccomplex_t IT_0133 = conj(IT_0068)*IT_0071;
    const ccomplex_t IT_0134 = (-6)*conj(IT_0032);
    const ccomplex_t IT_0135 = 0.166666666666667*IT_0016;
    const ccomplex_t IT_0136 = 6*s_34;
    const ccomplex_t IT_0137 = conj(IT_0068)*IT_0074;
    const ccomplex_t IT_0138 = (conj(IT_0013)*IT_0015 + IT_0016*conj(IT_0030))
      *IT_0033 + 6*IT_0068*((conj(IT_0013) + conj(IT_0032))*IT_0071 + (conj
      (IT_0030) + conj(IT_0072))*IT_0074 + 0.166666666666667*conj(IT_0068)*((-12
      )*IT_0014 + IT_0075)) + 6*IT_0108*(0.166666666666667*(12*IT_0014 + IT_0075
      )*conj(IT_0108) + IT_0109*(conj(IT_0030) + 0.166666666666667*IT_0110) +
       IT_0112*(conj(IT_0013) + IT_0113)) + -IT_0121*IT_0129 + IT_0124*(conj
      (IT_0013)*IT_0016 + IT_0014*conj(IT_0030) + IT_0130) + (-6)*IT_0032*(conj
      (IT_0108)*IT_0112 + -conj(IT_0072)*IT_0131 + (-0.166666666666667)*conj
      (IT_0032)*IT_0132 + -IT_0133) + 6*IT_0013*(conj(IT_0108)*IT_0112 + conj
      (IT_0030)*IT_0131 + 0.166666666666667*conj(IT_0013)*IT_0132 + IT_0133 +
       0.166666666666667*IT_0015*IT_0134 + IT_0110*IT_0135) + 6*IT_0030*
      (0.166666666666667*IT_0014*IT_0110 + IT_0130 + conj(IT_0013)*IT_0131 +
       IT_0134*IT_0135 + 0.166666666666667*conj(IT_0030)*IT_0136 + IT_0137) + 6
      *IT_0072*(conj(IT_0032)*IT_0131 + 0.166666666666667*conj(IT_0072)*IT_0136 
      + IT_0137);
    return create_ccomplex_return(IT_0138);
}

