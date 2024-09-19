#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_1_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_1_to_anti_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
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
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = (-0.333333333333333)*IT_0000;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_stau_1, 2);
    const ccomplex_t IT_0005 = cpow(s_12 + IT_0004 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = s_13*s_14;
    const ccomplex_t IT_0010 = pow(m_s, 2);
    const ccomplex_t IT_0011 = IT_0004*IT_0010;
    const ccomplex_t IT_0012 = s_34*IT_0004;
    const ccomplex_t IT_0013 = cos(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0014*IT_0015 +
       -IT_0013*IT_0016)*U_stau_00*conj(U_stau_00) + IT_0014*IT_0015*U_stau_10
      *conj(U_stau_10));
    const ccomplex_t IT_0018 = (-2)*IT_0017;
    const ccomplex_t IT_0019 = IT_0014*IT_0015;
    const ccomplex_t IT_0020 = e_em*IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0016;
    const ccomplex_t IT_0022 = e_em*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0020 + 3*IT_0022);
    const ccomplex_t IT_0024 = (-0.166666666666667)*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0004 + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0018*IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*IT_0014*IT_0015;
    const ccomplex_t IT_0031 = 0.333333333333333*IT_0030;
    const ccomplex_t IT_0032 = IT_0025*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0018*IT_0033;
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = (-6)*IT_0012;
    const ccomplex_t IT_0037 = (-6)*IT_0011;
    const ccomplex_t IT_0038 = 12*IT_0009;
    const ccomplex_t IT_0039 = IT_0036 + IT_0037 + IT_0038;
    const ccomplex_t IT_0040 = -s_14;
    const ccomplex_t IT_0041 = s_13 + IT_0040;
    const ccomplex_t IT_0042 = m_s*IT_0041;
    const ccomplex_t IT_0043 = -IT_0017;
    const ccomplex_t IT_0044 = m_s*IT_0043;
    const ccomplex_t IT_0045 = m_s*IT_0018;
    const ccomplex_t IT_0046 = -(IT_0027 + -IT_0033)*(IT_0044 + -IT_0045);
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = pow(m_W, -1);
    const ccomplex_t IT_0050 = cos(beta);
    const ccomplex_t IT_0051 = cpow(IT_0050, -1);
    const ccomplex_t IT_0052 = cos(alpha);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_s*e_em*IT_0016*IT_0049
      *IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = mu_h*m_tau*IT_0049;
    const ccomplex_t IT_0056 = sin(alpha);
    const ccomplex_t IT_0057 = sin(beta);
    const ccomplex_t IT_0058 = -IT_0056*IT_0057;
    const ccomplex_t IT_0059 = IT_0050*IT_0052;
    const ccomplex_t IT_0060 = IT_0058 + IT_0059;
    const ccomplex_t IT_0061 = cpow(IT_0013, -2);
    const ccomplex_t IT_0062 = pow(m_tau, 2);
    const ccomplex_t IT_0063 = (-2)*U_stau_00*(((-0.5)*m_W*IT_0060*(IT_0016 + 
      -IT_0015*IT_0061) + IT_0016*IT_0049*IT_0051*IT_0052*IT_0062)*conj
      (U_stau_00) + -1./2*IT_0016*IT_0051*IT_0055*IT_0056*conj(U_stau_10)) +
       U_stau_10*(IT_0016*IT_0051*IT_0055*IT_0056*conj(U_stau_00) + 2*(m_W
      *IT_0015*IT_0060*IT_0061 + -IT_0016*IT_0049*IT_0051*IT_0052*IT_0062)*conj
      (U_stau_10));
    const ccomplex_t IT_0064 = (0 + _Complex_I*-1)*(e_em*IT_0063 + (
      -1.4142135623731)*A_tau*IT_0052*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0054*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_12 + (-2)*IT_0004 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_s*e_em*IT_0016*IT_0049
      *IT_0051*IT_0056;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = IT_0050*IT_0056;
    const ccomplex_t IT_0073 = IT_0052*IT_0057;
    const ccomplex_t IT_0074 = IT_0072 + IT_0073;
    const ccomplex_t IT_0075 = 2*U_stau_00*((IT_0016*IT_0049*IT_0051*IT_0056
      *IT_0062 + (-0.5)*m_W*(IT_0016 + -IT_0015*IT_0061)*IT_0074)*conj(U_stau_00
      ) + 1./2*IT_0016*IT_0051*IT_0052*IT_0055*conj(U_stau_10)) + U_stau_10*
      (IT_0016*IT_0051*IT_0052*IT_0055*conj(U_stau_00) + 2*(IT_0016*IT_0049
      *IT_0051*IT_0056*IT_0062 + -m_W*IT_0015*IT_0061*IT_0074)*conj(U_stau_10));
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(e_em*IT_0075 +
       1.4142135623731*A_tau*IT_0056*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0071*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0004 + -reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = -IT_0069 + -IT_0081;
    const ccomplex_t IT_0083 = (IT_0027 + -IT_0033)*(IT_0044 + -IT_0045);
    const ccomplex_t IT_0084 = (-6)*IT_0048;
    const ccomplex_t IT_0085 = IT_0008*(24*IT_0009 + (-12)*IT_0011 + (-12)
      *IT_0012) + (IT_0029 + -conj(IT_0029) + IT_0035 + -conj(IT_0035))*IT_0039 
      + (-6)*IT_0042*(conj(IT_0048) + (-2)*IT_0082 + 2*conj(IT_0082) + -IT_0083 
      + conj(IT_0083) + 0.166666666666667*IT_0084);
    const ccomplex_t IT_0086 = m_s*e_em*IT_0016*IT_0049;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = e_em*mu_h*m_tau*IT_0049*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0089 = IT_0016*IT_0051*IT_0057;
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = e_em*mu_h*m_tau*IT_0049*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0092 = IT_0089*IT_0091;
    const ccomplex_t IT_0093 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0094 = IT_0050*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0097 = IT_0050*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (-0.5)*IT_0090 + 0.5*IT_0092 + 0.5*IT_0095 + (
      -0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0087*IT_0099;
    const ccomplex_t IT_0101 = pow(m_Z, 2);
    const ccomplex_t IT_0102 = cpow((-2)*s_12 + (-2)*IT_0004 + IT_0101 + 
      -reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0100*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = m_s*e_em*IT_0016*IT_0049*IT_0051*IT_0057;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = IT_0016*IT_0088;
    const ccomplex_t IT_0108 = IT_0016*IT_0091;
    const ccomplex_t IT_0109 = IT_0057*IT_0093;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = IT_0057*IT_0096;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (-0.5)*IT_0107 + 0.5*IT_0108 + (-0.5)*IT_0110 +
       0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0106*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_12 + (-2)*IT_0004 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = -IT_0104 + -IT_0117;
    const ccomplex_t IT_0119 = s_13 + s_14;
    const ccomplex_t IT_0120 = m_s*IT_0119;
    const ccomplex_t IT_0121 = s_34 + IT_0010;
    const ccomplex_t IT_0122 = (-6)*conj(IT_0048);
    const ccomplex_t IT_0123 = 12*s_34;
    const ccomplex_t IT_0124 = -IT_0010;
    const ccomplex_t IT_0125 = s_34 + IT_0124;
    const ccomplex_t IT_0126 = m_s*s_14;
    const ccomplex_t IT_0127 = conj(IT_0118)*IT_0121;
    const ccomplex_t IT_0128 = m_s*s_13;
    const ccomplex_t IT_0129 = IT_0036 + IT_0038;
    const ccomplex_t IT_0130 = IT_0042*conj(IT_0082);
    const ccomplex_t IT_0131 = (-6)*conj(IT_0029);
    const ccomplex_t IT_0132 = 0.166666666666667*IT_0126;
    const ccomplex_t IT_0133 = 6*s_34;
    const ccomplex_t IT_0134 = conj(IT_0082)*IT_0125;
    const ccomplex_t IT_0135 = -IT_0008*IT_0085 + (-6)*IT_0118*((conj(IT_0029)
       + -conj(IT_0035))*IT_0120 + -IT_0121*(conj(IT_0083) + 0.166666666666667
      *IT_0122) + (-0.166666666666667)*conj(IT_0118)*(12*IT_0010 + IT_0123)) + 6
      *IT_0082*((conj(IT_0029) + conj(IT_0035))*IT_0042 + 0.166666666666667*conj
      (IT_0082)*((-12)*IT_0010 + IT_0123) + (conj(IT_0048) + conj(IT_0083))
      *IT_0125) + IT_0084*(IT_0010*conj(IT_0083) + conj(IT_0035)*IT_0126 +
       IT_0127) + (-6)*IT_0029*(IT_0011*conj(IT_0035) + conj(IT_0118)*IT_0120 +
       conj(IT_0083)*IT_0126 + -conj(IT_0048)*IT_0128 + (-0.166666666666667)
      *conj(IT_0029)*IT_0129 + -IT_0130) + 6*IT_0035*(conj(IT_0118)*IT_0120 +
       conj(IT_0083)*IT_0128 + 0.166666666666667*conj(IT_0035)*IT_0129 + IT_0130
       + 0.166666666666667*IT_0011*IT_0131 + IT_0122*IT_0132) + 6*IT_0048*(conj
      (IT_0029)*IT_0128 + 0.166666666666667*conj(IT_0048)*IT_0133 + IT_0134) + 6
      *IT_0083*(0.166666666666667*IT_0010*IT_0122 + IT_0127 + conj(IT_0035)
      *IT_0128 + IT_0131*IT_0132 + 0.166666666666667*conj(IT_0083)*IT_0133 +
       IT_0134);
    return create_ccomplex_return(IT_0135);
}

