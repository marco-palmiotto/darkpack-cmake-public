#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_1_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_1_to_A0_Z(
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
      (U_stau_00) + -1./2*IT_0005*IT_0007*IT_0010*IT_0011*conj(U_stau_10)) +
       U_stau_10*(IT_0005*IT_0007*IT_0010*IT_0011*conj(U_stau_00) + 2*(m_W
      *IT_0006*IT_0014*IT_0015 + -IT_0002*IT_0007*IT_0009*IT_0011*IT_0016)*conj
      (U_stau_10));
    const ccomplex_t IT_0018 = (0 + _Complex_I*-1)*(e_em*IT_0017 + (
      -1.4142135623731)*A_tau*IT_0002*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_stau_1, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + (-2)*IT_0020 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0008*IT_0023;
    const ccomplex_t IT_0025 = e_em*IT_0001*(IT_0002*IT_0004 + IT_0003*IT_0005
      )*IT_0007;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = IT_0004*IT_0005;
    const ccomplex_t IT_0028 = IT_0002*IT_0003;
    const ccomplex_t IT_0029 = IT_0027 + IT_0028;
    const ccomplex_t IT_0030 = 2*U_stau_00*((IT_0005*IT_0007*IT_0009*IT_0011
      *IT_0016 + (-0.5)*m_W*(IT_0007 + -IT_0006*IT_0015)*IT_0029)*conj(U_stau_00
      ) + 1./2*IT_0002*IT_0007*IT_0010*IT_0011*conj(U_stau_10)) + U_stau_10*
      (IT_0002*IT_0007*IT_0010*IT_0011*conj(U_stau_00) + 2*(IT_0005*IT_0007
      *IT_0009*IT_0011*IT_0016 + -m_W*IT_0006*IT_0015*IT_0029)*conj(U_stau_10));
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(e_em*IT_0030 +
       1.4142135623731*A_tau*IT_0005*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0020 + -reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0026*IT_0035;
    const ccomplex_t IT_0037 = -IT_0024 + -IT_0036;
    const ccomplex_t IT_0038 = pow(m_Z, -2);
    const ccomplex_t IT_0039 = pow(s_14, 2);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0006 +
       -IT_0000*IT_0007)*U_stau_00*conj(U_stau_01) + IT_0001*IT_0006*U_stau_10
      *conj(U_stau_11));
    const ccomplex_t IT_0042 = (-2)*IT_0041;
    const ccomplex_t IT_0043 = e_em*mu_h*m_tau*IT_0009*U_stau_01*conj
      (U_stau_10);
    const ccomplex_t IT_0044 = IT_0007*IT_0043;
    const ccomplex_t IT_0045 = e_em*mu_h*m_tau*IT_0009*conj(U_stau_00)
      *U_stau_11;
    const ccomplex_t IT_0046 = IT_0007*IT_0045;
    const ccomplex_t IT_0047 = A_tau*U_stau_01*conj(U_stau_10);
    const ccomplex_t IT_0048 = IT_0003*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = A_tau*conj(U_stau_00)*U_stau_11;
    const ccomplex_t IT_0051 = IT_0003*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = 0.5*IT_0044 + (-0.5)*IT_0046 + 0.5*IT_0049 + (
      -0.5)*IT_0052;
    const ccomplex_t IT_0054 = pow(m_A0, 2);
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0020 + IT_0054 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0042*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0006 +
       -IT_0000*IT_0007)*U_stau_00*conj(U_stau_00) + IT_0001*IT_0006*U_stau_10
      *conj(U_stau_10));
    const ccomplex_t IT_0060 = (-2)*IT_0059;
    const ccomplex_t IT_0061 = e_em*mu_h*m_tau*IT_0009*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0062 = IT_0007*IT_0061;
    const ccomplex_t IT_0063 = e_em*mu_h*m_tau*IT_0009*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0064 = IT_0007*IT_0063;
    const ccomplex_t IT_0065 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0066 = IT_0003*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0069 = IT_0003*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (-0.5)*IT_0062 + 0.5*IT_0064 + (-0.5)*IT_0067 +
       0.5*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_23 + IT_0054 + reg_prop + (0 +
       _Complex_I*1)*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0060*IT_0074;
    const ccomplex_t IT_0076 = -IT_0058 + -IT_0075;
    const ccomplex_t IT_0077 = s_14*s_34*IT_0038;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = s_13 + IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0006 +
       -IT_0000*IT_0007)*conj(U_stau_00)*U_stau_01 + IT_0001*IT_0006*conj
      (U_stau_10)*U_stau_11);
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = e_em*mu_h*m_tau*IT_0009*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0083 = IT_0007*IT_0082;
    const ccomplex_t IT_0084 = e_em*mu_h*m_tau*IT_0009*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0085 = IT_0007*IT_0084;
    const ccomplex_t IT_0086 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0087 = IT_0003*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0090 = IT_0003*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (-0.5)*IT_0083 + 0.5*IT_0085 + (-0.5)*IT_0088 +
       0.5*IT_0091;
    const ccomplex_t IT_0093 = cpow((-2)*s_13 + IT_0020 + IT_0054 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0081*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_13 + IT_0054 + reg_prop + (0 +
       _Complex_I*1)*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0098 = IT_0071*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = IT_0060*IT_0099;
    const ccomplex_t IT_0101 = -IT_0024 + -IT_0036 + IT_0058 + IT_0075 +
       IT_0096 + IT_0100;
    const ccomplex_t IT_0102 = s_14*s_24*IT_0038;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = s_12 + IT_0103;
    const ccomplex_t IT_0105 = -IT_0037;
    const ccomplex_t IT_0106 = -conj(IT_0037);
    const ccomplex_t IT_0107 = pow(s_34, 2);
    const ccomplex_t IT_0108 = IT_0038*IT_0107;
    const ccomplex_t IT_0109 = s_24*s_34*IT_0038;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = s_23 + IT_0110;
    const ccomplex_t IT_0112 = pow(s_24, 2);
    const ccomplex_t IT_0113 = IT_0038*IT_0112;
    const ccomplex_t IT_0114 = -IT_0101;
    const ccomplex_t IT_0115 = (conj(IT_0037)*(IT_0020 + -IT_0040) + conj
      (IT_0076)*IT_0079 + conj(IT_0101)*IT_0104)*IT_0105 + IT_0101*IT_0104
      *IT_0106 + IT_0076*(IT_0079*IT_0106 + -conj(IT_0076)*(IT_0054 + -IT_0108) 
      + -conj(IT_0101)*IT_0111) + (conj(IT_0076)*IT_0111 + conj(IT_0101)*
      (IT_0020 + -IT_0113))*IT_0114;
    return create_ccomplex_return(IT_0115);
}

