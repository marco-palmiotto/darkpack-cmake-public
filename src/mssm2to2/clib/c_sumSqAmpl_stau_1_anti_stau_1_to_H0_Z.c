#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_1_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_1_to_H0_Z(
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
       -IT_0000*IT_0003)*U_stau_00*conj(U_stau_00) + IT_0001*IT_0002*U_stau_10
      *conj(U_stau_10));
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
    const ccomplex_t IT_0025 = pow(m_stau_1, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0025 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0005*IT_0028;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = pow(m_Z, -2);
    const ccomplex_t IT_0032 = pow(s_14, 2);
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = (-2)*IT_0004;
    const ccomplex_t IT_0035 = pow(m_W, -1);
    const ccomplex_t IT_0036 = mu_h*m_tau*IT_0035;
    const ccomplex_t IT_0037 = cpow(IT_0006, -1);
    const ccomplex_t IT_0038 = -IT_0010*IT_0011;
    const ccomplex_t IT_0039 = IT_0006*IT_0007;
    const ccomplex_t IT_0040 = IT_0038 + IT_0039;
    const ccomplex_t IT_0041 = pow(m_tau, 2);
    const ccomplex_t IT_0042 = (-2)*U_stau_00*(((-0.5)*m_W*(IT_0003 + -IT_0002
      *IT_0016)*IT_0040 + IT_0003*IT_0007*IT_0035*IT_0037*IT_0041)*conj
      (U_stau_00) + -1./2*IT_0003*IT_0011*IT_0036*IT_0037*conj(U_stau_10)) +
       U_stau_10*(IT_0003*IT_0011*IT_0036*IT_0037*conj(U_stau_00) + 2*(m_W
      *IT_0002*IT_0016*IT_0040 + -IT_0003*IT_0007*IT_0035*IT_0037*IT_0041)*conj
      (U_stau_10));
    const ccomplex_t IT_0043 = (0 + _Complex_I*-1)*(e_em*IT_0042 + (
      -1.4142135623731)*A_tau*IT_0007*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = pow(m_H0, 2);
    const ccomplex_t IT_0046 = cpow((-2)*s_23 + IT_0045 + reg_prop + (0 +
       _Complex_I*1)*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0034*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*U_stau_00*conj(U_stau_01) + IT_0001*IT_0002*U_stau_10
      *conj(U_stau_11));
    const ccomplex_t IT_0051 = (-2)*IT_0050;
    const ccomplex_t IT_0052 = (-2)*conj(U_stau_00)*(((-0.5)*m_W*(IT_0003 + 
      -IT_0002*IT_0016)*IT_0040 + IT_0003*IT_0007*IT_0035*IT_0037*IT_0041)
      *U_stau_01 + -1./2*IT_0003*IT_0011*IT_0036*IT_0037*U_stau_11) + conj
      (U_stau_10)*(IT_0003*IT_0011*IT_0036*IT_0037*U_stau_01 + 2*(m_W*IT_0002
      *IT_0016*IT_0040 + -IT_0003*IT_0007*IT_0035*IT_0037*IT_0041)*U_stau_11);
    const ccomplex_t IT_0053 = (0 + _Complex_I*-1)*(e_em*IT_0052 + (
      -1.4142135623731)*A_tau*IT_0007*(U_stau_01*conj(U_stau_10) + conj
      (U_stau_00)*U_stau_11));
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0025 + IT_0045 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0051*IT_0057;
    const ccomplex_t IT_0059 = e_em*IT_0001*IT_0003*(IT_0007*IT_0010 + 
      -IT_0006*IT_0011);
    const ccomplex_t IT_0060 = e_em*mu_h*m_tau*IT_0035*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0061 = IT_0003*IT_0060;
    const ccomplex_t IT_0062 = e_em*mu_h*m_tau*IT_0035*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0065 = IT_0010*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0068 = IT_0010*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (-0.5)*IT_0061 + 0.5*IT_0063 + (-0.5)*IT_0066 +
       0.5*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_12 + (-2)*IT_0025 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0059*IT_0073;
    const ccomplex_t IT_0075 = e_em*IT_0001*IT_0003*(IT_0006*IT_0007 + IT_0010
      *IT_0011);
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0003*IT_0010*IT_0037;
    const ccomplex_t IT_0078 = IT_0060*IT_0077;
    const ccomplex_t IT_0079 = IT_0062*IT_0077;
    const ccomplex_t IT_0080 = IT_0006*IT_0064;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = IT_0006*IT_0067;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (-0.5)*IT_0078 + 0.5*IT_0079 + 0.5*IT_0081 + (
      -0.5)*IT_0083;
    const ccomplex_t IT_0085 = pow(m_Z, 2);
    const ccomplex_t IT_0086 = cpow((-2)*s_12 + (-2)*IT_0025 + IT_0085 + 
      -reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0076*IT_0088;
    const ccomplex_t IT_0090 = -IT_0049 + -IT_0058 + IT_0074 + IT_0089;
    const ccomplex_t IT_0091 = s_14*s_34*IT_0031;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = s_13 + IT_0092;
    const ccomplex_t IT_0094 = -IT_0034;
    const ccomplex_t IT_0095 = IT_0005 + IT_0094;
    const ccomplex_t IT_0096 = IT_0028*IT_0095;
    const ccomplex_t IT_0097 = (-2)*U_stau_00*(((-0.5)*m_W*(IT_0003 + -IT_0002
      *IT_0016)*IT_0040 + IT_0003*IT_0007*IT_0035*IT_0037*IT_0041)*conj
      (U_stau_00) + -1./2*IT_0003*IT_0011*IT_0036*IT_0037*conj(U_stau_10)) +
       U_stau_10*(IT_0003*IT_0011*IT_0036*IT_0037*conj(U_stau_00) + 2*(m_W
      *IT_0002*IT_0016*IT_0040 + -IT_0003*IT_0007*IT_0035*IT_0037*IT_0041)*conj
      (U_stau_10));
    const ccomplex_t IT_0098 = (0 + _Complex_I*-1)*(e_em*IT_0097 + (
      -1.4142135623731)*A_tau*IT_0007*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = cpow((-2)*s_13 + IT_0045 + reg_prop + (0 +
       _Complex_I*1)*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0101 = IT_0099*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0034*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*conj(U_stau_00)*U_stau_01 + IT_0001*IT_0002*conj
      (U_stau_10)*U_stau_11);
    const ccomplex_t IT_0105 = (-2)*IT_0104;
    const ccomplex_t IT_0106 = (-2)*U_stau_00*(((-0.5)*m_W*(IT_0003 + -IT_0002
      *IT_0016)*IT_0040 + IT_0003*IT_0007*IT_0035*IT_0037*IT_0041)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0011*IT_0036*IT_0037*conj(U_stau_11)) +
       U_stau_10*(IT_0003*IT_0011*IT_0036*IT_0037*conj(U_stau_01) + 2*(m_W
      *IT_0002*IT_0016*IT_0040 + -IT_0003*IT_0007*IT_0035*IT_0037*IT_0041)*conj
      (U_stau_11));
    const ccomplex_t IT_0107 = (0 + _Complex_I*-1)*(e_em*IT_0106 + (
      -1.4142135623731)*A_tau*IT_0007*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_13 + IT_0025 + IT_0045 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0110 = IT_0108*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0105*IT_0111;
    const ccomplex_t IT_0113 = IT_0049 + IT_0058 + -IT_0096 + IT_0103 + IT_0112;
    const ccomplex_t IT_0114 = s_14*s_24*IT_0031;
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = s_12 + IT_0115;
    const ccomplex_t IT_0117 = -IT_0030;
    const ccomplex_t IT_0118 = pow(s_34, 2);
    const ccomplex_t IT_0119 = IT_0031*IT_0118;
    const ccomplex_t IT_0120 = s_24*s_34*IT_0031;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = s_23 + IT_0121;
    const ccomplex_t IT_0123 = -conj(IT_0030);
    const ccomplex_t IT_0124 = pow(s_24, 2);
    const ccomplex_t IT_0125 = IT_0031*IT_0124;
    const ccomplex_t IT_0126 = (conj(IT_0030)*(IT_0025 + -IT_0033) + conj
      (IT_0090)*IT_0093 + conj(IT_0113)*IT_0116)*IT_0117 + -IT_0090*(conj
      (IT_0090)*(IT_0045 + -IT_0119) + conj(IT_0113)*IT_0122 + -IT_0093*IT_0123)
       + -IT_0113*(conj(IT_0090)*IT_0122 + -IT_0116*IT_0123 + conj(IT_0113)*
      (IT_0025 + -IT_0125));
    return create_ccomplex_return(IT_0126);
}

