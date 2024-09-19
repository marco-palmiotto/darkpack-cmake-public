#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_stau_2_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_stau_2_to_H0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
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
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
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
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*U_stau_01*conj(U_stau_01) + IT_0001*IT_0002*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cos(alpha);
    const ccomplex_t IT_0008 = cpow(IT_0000, -4);
    const ccomplex_t IT_0009 = IT_0002*IT_0006*IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = IT_0002*IT_0008*IT_0010*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0002, -3);
    const ccomplex_t IT_0014 = IT_0006*IT_0007*IT_0013;
    const ccomplex_t IT_0015 = IT_0010*IT_0011*IT_0013;
    const ccomplex_t IT_0016 = cpow(IT_0000, -2);
    const ccomplex_t IT_0017 = IT_0003*IT_0006*IT_0007*IT_0016;
    const ccomplex_t IT_0018 = IT_0003*IT_0010*IT_0011*IT_0016;
    const ccomplex_t IT_0019 = cpow(IT_0002, -2);
    const ccomplex_t IT_0020 = cpow(IT_0016 + IT_0019, -1);
    const ccomplex_t IT_0021 = (IT_0009 + IT_0012 + IT_0014 + IT_0015 + 2
      *IT_0017 + 2*IT_0018)*IT_0020;
    const ccomplex_t IT_0022 = m_W*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = pow(m_stau_2, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0025 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0005*IT_0028;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = pow(m_Z, -2);
    const ccomplex_t IT_0032 = pow(s_14, 2);
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = e_em*IT_0001*IT_0003*(IT_0006*IT_0007 + IT_0010
      *IT_0011);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = pow(m_W, -1);
    const ccomplex_t IT_0037 = e_em*mu_h*m_tau*IT_0036*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0038 = cpow(IT_0006, -1);
    const ccomplex_t IT_0039 = IT_0003*IT_0010*IT_0038;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = e_em*mu_h*m_tau*IT_0036*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0044 = IT_0006*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0047 = IT_0006*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (-0.5)*IT_0040 + 0.5*IT_0042 + 0.5*IT_0045 + (
      -0.5)*IT_0048;
    const ccomplex_t IT_0050 = pow(m_Z, 2);
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + (-2)*IT_0025 + IT_0050 + 
      -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0035*IT_0053;
    const ccomplex_t IT_0055 = (-2)*IT_0004;
    const ccomplex_t IT_0056 = mu_h*m_tau*IT_0036;
    const ccomplex_t IT_0057 = -IT_0010*IT_0011;
    const ccomplex_t IT_0058 = IT_0006*IT_0007;
    const ccomplex_t IT_0059 = IT_0057 + IT_0058;
    const ccomplex_t IT_0060 = pow(m_tau, 2);
    const ccomplex_t IT_0061 = (-2)*U_stau_01*(((-0.5)*m_W*(IT_0003 + -IT_0002
      *IT_0016)*IT_0059 + IT_0003*IT_0007*IT_0036*IT_0038*IT_0060)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0011*IT_0038*IT_0056*conj(U_stau_11)) +
       U_stau_11*(IT_0003*IT_0011*IT_0038*IT_0056*conj(U_stau_01) + 2*(m_W
      *IT_0002*IT_0016*IT_0059 + -IT_0003*IT_0007*IT_0036*IT_0038*IT_0060)*conj
      (U_stau_11));
    const ccomplex_t IT_0062 = (0 + _Complex_I*-1)*(e_em*IT_0061 + (
      -1.4142135623731)*A_tau*IT_0007*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = pow(m_H0, 2);
    const ccomplex_t IT_0065 = cpow((-2)*s_23 + IT_0064 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0055*IT_0067;
    const ccomplex_t IT_0069 = e_em*IT_0001*IT_0003*(IT_0007*IT_0010 + 
      -IT_0006*IT_0011);
    const ccomplex_t IT_0070 = IT_0003*IT_0037;
    const ccomplex_t IT_0071 = IT_0003*IT_0041;
    const ccomplex_t IT_0072 = IT_0010*IT_0043;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = IT_0010*IT_0046;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = (-0.5)*IT_0070 + 0.5*IT_0071 + (-0.5)*IT_0073 +
       0.5*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_12 + (-2)*IT_0025 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0069*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*conj(U_stau_00)*U_stau_01 + IT_0001*IT_0002*conj
      (U_stau_10)*U_stau_11);
    const ccomplex_t IT_0082 = (-2)*IT_0081;
    const ccomplex_t IT_0083 = (-2)*U_stau_00*(((-0.5)*m_W*(IT_0003 + -IT_0002
      *IT_0016)*IT_0059 + IT_0003*IT_0007*IT_0036*IT_0038*IT_0060)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0011*IT_0038*IT_0056*conj(U_stau_11)) +
       U_stau_10*(IT_0003*IT_0011*IT_0038*IT_0056*conj(U_stau_01) + 2*(m_W
      *IT_0002*IT_0016*IT_0059 + -IT_0003*IT_0007*IT_0036*IT_0038*IT_0060)*conj
      (U_stau_11));
    const ccomplex_t IT_0084 = (0 + _Complex_I*-1)*(e_em*IT_0083 + (
      -1.4142135623731)*A_tau*IT_0007*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0025 + IT_0064 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0082*IT_0088;
    const ccomplex_t IT_0090 = IT_0054 + -IT_0068 + IT_0080 + -IT_0089;
    const ccomplex_t IT_0091 = s_14*s_34*IT_0031;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = s_13 + IT_0092;
    const ccomplex_t IT_0094 = -IT_0055;
    const ccomplex_t IT_0095 = IT_0005 + IT_0094;
    const ccomplex_t IT_0096 = IT_0028*IT_0095;
    const ccomplex_t IT_0097 = (-2)*U_stau_01*(((-0.5)*m_W*(IT_0003 + -IT_0002
      *IT_0016)*IT_0059 + IT_0003*IT_0007*IT_0036*IT_0038*IT_0060)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0011*IT_0038*IT_0056*conj(U_stau_11)) +
       U_stau_11*(IT_0003*IT_0011*IT_0038*IT_0056*conj(U_stau_01) + 2*(m_W
      *IT_0002*IT_0016*IT_0059 + -IT_0003*IT_0007*IT_0036*IT_0038*IT_0060)*conj
      (U_stau_11));
    const ccomplex_t IT_0098 = (0 + _Complex_I*-1)*(e_em*IT_0097 + (
      -1.4142135623731)*A_tau*IT_0007*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = cpow((-2)*s_13 + IT_0064 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0101 = IT_0099*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0055*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*U_stau_00*conj(U_stau_01) + IT_0001*IT_0002*U_stau_10
      *conj(U_stau_11));
    const ccomplex_t IT_0105 = (-2)*IT_0104;
    const ccomplex_t IT_0106 = (-2)*conj(U_stau_00)*(((-0.5)*m_W*(IT_0003 + 
      -IT_0002*IT_0016)*IT_0059 + IT_0003*IT_0007*IT_0036*IT_0038*IT_0060)
      *U_stau_01 + -1./2*IT_0003*IT_0011*IT_0038*IT_0056*U_stau_11) + conj
      (U_stau_10)*(IT_0003*IT_0011*IT_0038*IT_0056*U_stau_01 + 2*(m_W*IT_0002
      *IT_0016*IT_0059 + -IT_0003*IT_0007*IT_0036*IT_0038*IT_0060)*U_stau_11);
    const ccomplex_t IT_0107 = (0 + _Complex_I*-1)*(e_em*IT_0106 + (
      -1.4142135623731)*A_tau*IT_0007*(U_stau_01*conj(U_stau_10) + conj
      (U_stau_00)*U_stau_11));
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_13 + IT_0025 + IT_0064 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0110 = IT_0108*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0105*IT_0111;
    const ccomplex_t IT_0113 = (-2)*U_stau_00*(((-0.5)*m_W*(IT_0003 + -IT_0002
      *IT_0016)*IT_0059 + IT_0003*IT_0007*IT_0036*IT_0038*IT_0060)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0011*IT_0038*IT_0056*conj(U_stau_11)) +
       U_stau_10*(IT_0003*IT_0011*IT_0038*IT_0056*conj(U_stau_01) + 2*(m_W
      *IT_0002*IT_0016*IT_0059 + -IT_0003*IT_0007*IT_0036*IT_0038*IT_0060)*conj
      (U_stau_11));
    const ccomplex_t IT_0114 = (0 + _Complex_I*-1)*(e_em*IT_0113 + (
      -1.4142135623731)*A_tau*IT_0007*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = IT_0086*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = IT_0082*IT_0117;
    const ccomplex_t IT_0119 = IT_0068 + -IT_0096 + IT_0103 + IT_0112 + IT_0118;
    const ccomplex_t IT_0120 = s_14*s_24*IT_0031;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = s_12 + IT_0121;
    const ccomplex_t IT_0123 = -IT_0030;
    const ccomplex_t IT_0124 = -conj(IT_0030);
    const ccomplex_t IT_0125 = pow(s_34, 2);
    const ccomplex_t IT_0126 = IT_0031*IT_0125;
    const ccomplex_t IT_0127 = s_24*s_34*IT_0031;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = s_23 + IT_0128;
    const ccomplex_t IT_0130 = pow(s_24, 2);
    const ccomplex_t IT_0131 = IT_0031*IT_0130;
    const ccomplex_t IT_0132 = -IT_0119;
    const ccomplex_t IT_0133 = (conj(IT_0030)*(IT_0025 + -IT_0033) + conj
      (IT_0090)*IT_0093 + conj(IT_0119)*IT_0122)*IT_0123 + IT_0119*IT_0122
      *IT_0124 + IT_0090*(IT_0093*IT_0124 + -conj(IT_0090)*(IT_0064 + -IT_0126) 
      + -conj(IT_0119)*IT_0129) + (conj(IT_0090)*IT_0129 + conj(IT_0119)*
      (IT_0025 + -IT_0131))*IT_0132;
    return create_ccomplex_return(IT_0133);
}

