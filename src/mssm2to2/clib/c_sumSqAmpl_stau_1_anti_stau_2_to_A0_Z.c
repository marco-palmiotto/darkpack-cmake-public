#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_2_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_2_to_A0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
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
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = e_em*IT_0001*(IT_0002*IT_0003 + -IT_0004
      *IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = mu_h*m_tau*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0004, -1);
    const ccomplex_t IT_0012 = -IT_0003*IT_0005;
    const ccomplex_t IT_0013 = IT_0002*IT_0004;
    const ccomplex_t IT_0014 = IT_0012 + IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0000, -2);
    const ccomplex_t IT_0016 = pow(m_tau, 2);
    const ccomplex_t IT_0017 = (-2)*U_stau_00*(((-0.5)*m_W*IT_0014*(IT_0007 + 
      -IT_0006*IT_0015) + IT_0002*IT_0007*IT_0009*IT_0011*IT_0016)*conj
      (U_stau_01) + -1./2*IT_0005*IT_0007*IT_0010*IT_0011*conj(U_stau_11)) +
       U_stau_10*(IT_0005*IT_0007*IT_0010*IT_0011*conj(U_stau_01) + 2*(m_W
      *IT_0006*IT_0014*IT_0015 + -IT_0002*IT_0007*IT_0009*IT_0011*IT_0016)*conj
      (U_stau_11));
    const ccomplex_t IT_0018 = (0 + _Complex_I*-1)*(e_em*IT_0017 + (
      -1.4142135623731)*A_tau*IT_0002*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_stau_1, 2);
    const ccomplex_t IT_0021 = pow(m_stau_2, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + -IT_0020 + -IT_0021 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0008*IT_0024;
    const ccomplex_t IT_0026 = e_em*IT_0001*(IT_0002*IT_0004 + IT_0003*IT_0005
      )*IT_0007;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0004*IT_0005;
    const ccomplex_t IT_0029 = IT_0002*IT_0003;
    const ccomplex_t IT_0030 = IT_0028 + IT_0029;
    const ccomplex_t IT_0031 = 2*U_stau_00*((IT_0005*IT_0007*IT_0009*IT_0011
      *IT_0016 + (-0.5)*m_W*(IT_0007 + -IT_0006*IT_0015)*IT_0030)*conj(U_stau_01
      ) + 1./2*IT_0002*IT_0007*IT_0010*IT_0011*conj(U_stau_11)) + U_stau_10*
      (IT_0002*IT_0007*IT_0010*IT_0011*conj(U_stau_01) + 2*(IT_0005*IT_0007
      *IT_0009*IT_0011*IT_0016 + -m_W*IT_0006*IT_0015*IT_0030)*conj(U_stau_11));
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(e_em*IT_0031 +
       1.4142135623731*A_tau*IT_0005*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0020 + -IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = IT_0027*IT_0036;
    const ccomplex_t IT_0038 = -IT_0025 + -IT_0037;
    const ccomplex_t IT_0039 = pow(m_Z, -2);
    const ccomplex_t IT_0040 = pow(s_14, 2);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0006 +
       -IT_0000*IT_0007)*U_stau_00*conj(U_stau_00) + IT_0001*IT_0006*U_stau_10
      *conj(U_stau_10));
    const ccomplex_t IT_0043 = (-2)*IT_0042;
    const ccomplex_t IT_0044 = e_em*mu_h*m_tau*IT_0009*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0045 = IT_0007*IT_0044;
    const ccomplex_t IT_0046 = e_em*mu_h*m_tau*IT_0009*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0047 = IT_0007*IT_0046;
    const ccomplex_t IT_0048 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0049 = IT_0003*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0052 = IT_0003*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (-0.5)*IT_0045 + 0.5*IT_0047 + (-0.5)*IT_0050 +
       0.5*IT_0053;
    const ccomplex_t IT_0055 = pow(m_A0, 2);
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0021 + IT_0055 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0043*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0006 +
       -IT_0000*IT_0007)*U_stau_00*conj(U_stau_01) + IT_0001*IT_0006*U_stau_10
      *conj(U_stau_11));
    const ccomplex_t IT_0061 = (-2)*IT_0060;
    const ccomplex_t IT_0062 = e_em*mu_h*m_tau*IT_0009*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0063 = IT_0007*IT_0062;
    const ccomplex_t IT_0064 = e_em*mu_h*m_tau*IT_0009*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0065 = IT_0007*IT_0064;
    const ccomplex_t IT_0066 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0067 = IT_0003*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0070 = IT_0003*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (-0.5)*IT_0063 + 0.5*IT_0065 + (-0.5)*IT_0068 +
       0.5*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_23 + IT_0055 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0061*IT_0075;
    const ccomplex_t IT_0077 = -IT_0059 + -IT_0076;
    const ccomplex_t IT_0078 = s_14*s_34*IT_0039;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = s_13 + IT_0079;
    const ccomplex_t IT_0081 = e_em*mu_h*m_tau*IT_0009*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0082 = IT_0007*IT_0081;
    const ccomplex_t IT_0083 = e_em*mu_h*m_tau*IT_0009*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0084 = IT_0007*IT_0083;
    const ccomplex_t IT_0085 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0086 = IT_0003*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0089 = IT_0003*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0082 + 0.5*IT_0084 + (-0.5)*IT_0087 +
       0.5*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_13 + IT_0055 + reg_prop + (0 +
       _Complex_I*1)*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = IT_0061*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0006 +
       -IT_0000*IT_0007)*U_stau_01*conj(U_stau_01) + IT_0001*IT_0006*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0097 = (-2)*IT_0096;
    const ccomplex_t IT_0098 = cpow((-2)*s_13 + IT_0020 + IT_0055 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0099 = IT_0054*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0097*IT_0100;
    const ccomplex_t IT_0102 = 2*U_stau_00*((IT_0005*IT_0007*IT_0009*IT_0011
      *IT_0016 + (-0.5)*m_W*(IT_0007 + -IT_0006*IT_0015)*IT_0030)*conj(U_stau_01
      ) + 1./2*IT_0002*IT_0007*IT_0010*IT_0011*conj(U_stau_11)) + U_stau_10*
      (IT_0002*IT_0007*IT_0010*IT_0011*conj(U_stau_01) + 2*(IT_0005*IT_0007
      *IT_0009*IT_0011*IT_0016 + -m_W*IT_0006*IT_0015*IT_0030)*conj(U_stau_11));
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(e_em*IT_0102 +
       1.4142135623731*A_tau*IT_0005*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0034*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0027*IT_0106;
    const ccomplex_t IT_0108 = (-2)*U_stau_00*(((-0.5)*m_W*IT_0014*(IT_0007 + 
      -IT_0006*IT_0015) + IT_0002*IT_0007*IT_0009*IT_0011*IT_0016)*conj
      (U_stau_01) + -1./2*IT_0005*IT_0007*IT_0010*IT_0011*conj(U_stau_11)) +
       U_stau_10*(IT_0005*IT_0007*IT_0010*IT_0011*conj(U_stau_01) + 2*(m_W
      *IT_0006*IT_0014*IT_0015 + -IT_0002*IT_0007*IT_0009*IT_0011*IT_0016)*conj
      (U_stau_11));
    const ccomplex_t IT_0109 = (0 + _Complex_I*-1)*(e_em*IT_0108 + (
      -1.4142135623731)*A_tau*IT_0002*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0110 = (-0.5)*IT_0109;
    const ccomplex_t IT_0111 = IT_0022*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = IT_0008*IT_0112;
    const ccomplex_t IT_0114 = IT_0059 + IT_0076 + IT_0095 + IT_0101 + 
      -IT_0107 + -IT_0113;
    const ccomplex_t IT_0115 = s_14*s_24*IT_0039;
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = s_12 + IT_0116;
    const ccomplex_t IT_0118 = -IT_0038;
    const ccomplex_t IT_0119 = -conj(IT_0038);
    const ccomplex_t IT_0120 = pow(s_34, 2);
    const ccomplex_t IT_0121 = IT_0039*IT_0120;
    const ccomplex_t IT_0122 = s_24*s_34*IT_0039;
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = s_23 + IT_0123;
    const ccomplex_t IT_0125 = pow(s_24, 2);
    const ccomplex_t IT_0126 = IT_0039*IT_0125;
    const ccomplex_t IT_0127 = -IT_0114;
    const ccomplex_t IT_0128 = (conj(IT_0038)*(IT_0020 + -IT_0041) + conj
      (IT_0077)*IT_0080 + conj(IT_0114)*IT_0117)*IT_0118 + IT_0114*IT_0117
      *IT_0119 + IT_0077*(IT_0080*IT_0119 + -conj(IT_0077)*(IT_0055 + -IT_0121) 
      + -conj(IT_0114)*IT_0124) + (conj(IT_0077)*IT_0124 + conj(IT_0114)*
      (IT_0021 + -IT_0126))*IT_0127;
    return create_ccomplex_return(IT_0128);
}

