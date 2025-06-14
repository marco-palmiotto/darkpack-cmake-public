#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_2_to_anti_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_2_to_anti_mu_mu(
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
    const creal_t m_mu = param->m_mu;
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
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*m_mu*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = mu_h*m_tau*IT_0000;
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = -IT_0009*IT_0010;
    const ccomplex_t IT_0012 = IT_0001*IT_0003;
    const ccomplex_t IT_0013 = IT_0011 + IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -2);
    const ccomplex_t IT_0016 = pow(m_tau, 2);
    const ccomplex_t IT_0017 = (-2)*U_stau_00*(((-0.5)*m_W*IT_0013*(IT_0005 + 
      -IT_0004*IT_0015) + IT_0000*IT_0002*IT_0003*IT_0005*IT_0016)*conj
      (U_stau_01) + -1./2*IT_0002*IT_0005*IT_0008*IT_0009*conj(U_stau_11)) +
       U_stau_10*(IT_0002*IT_0005*IT_0008*IT_0009*conj(U_stau_01) + 2*(m_W
      *IT_0004*IT_0013*IT_0015 + -IT_0000*IT_0002*IT_0003*IT_0005*IT_0016)*conj
      (U_stau_11));
    const ccomplex_t IT_0018 = (0 + _Complex_I*-1)*(e_em*IT_0017 + (
      -1.4142135623731)*A_tau*IT_0003*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = pow(m_stau_1, 2);
    const ccomplex_t IT_0022 = pow(m_stau_2, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + -IT_0021 + -IT_0022 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*e_em*m_mu*IT_0000*IT_0002
      *IT_0005*IT_0009;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0001*IT_0009;
    const ccomplex_t IT_0029 = IT_0003*IT_0010;
    const ccomplex_t IT_0030 = IT_0028 + IT_0029;
    const ccomplex_t IT_0031 = 2*U_stau_00*((IT_0000*IT_0002*IT_0005*IT_0009
      *IT_0016 + (-0.5)*m_W*(IT_0005 + -IT_0004*IT_0015)*IT_0030)*conj(U_stau_01
      ) + 1./2*IT_0002*IT_0003*IT_0005*IT_0008*conj(U_stau_11)) + U_stau_10*
      (IT_0002*IT_0003*IT_0005*IT_0008*conj(U_stau_01) + 2*(IT_0000*IT_0002
      *IT_0005*IT_0009*IT_0016 + -m_W*IT_0004*IT_0015*IT_0030)*conj(U_stau_11));
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(e_em*IT_0031 +
       1.4142135623731*A_tau*IT_0009*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0021 + -IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -IT_0025 + -IT_0037;
    const ccomplex_t IT_0039 = cpow(IT_0014, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0005
      *IT_0014 + -IT_0004*IT_0039)*U_stau_00*conj(U_stau_01) + IT_0004*IT_0039
      *U_stau_10*conj(U_stau_11));
    const ccomplex_t IT_0041 = (-2)*IT_0040;
    const ccomplex_t IT_0042 = IT_0004*IT_0039;
    const ccomplex_t IT_0043 = e_em*IT_0042;
    const ccomplex_t IT_0044 = IT_0005*IT_0014;
    const ccomplex_t IT_0045 = e_em*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0043 + -IT_0045);
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0021 + -IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0041*IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0039;
    const ccomplex_t IT_0054 = IT_0048*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0041*IT_0055;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = -s_14;
    const ccomplex_t IT_0059 = s_13 + IT_0058;
    const ccomplex_t IT_0060 = m_mu*IT_0059;
    const ccomplex_t IT_0061 = -IT_0040;
    const ccomplex_t IT_0062 = m_mu*IT_0061;
    const ccomplex_t IT_0063 = m_mu*IT_0041;
    const ccomplex_t IT_0064 = -(IT_0050 + -IT_0055)*(IT_0062 + -IT_0063);
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = pow(m_mu, 2);
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = s_34 + IT_0068;
    const ccomplex_t IT_0070 = m_mu*s_13;
    const ccomplex_t IT_0071 = e_em*m_mu*IT_0000*IT_0002*IT_0005*IT_0010;
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = e_em*mu_h*m_tau*IT_0000*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0074 = IT_0005*IT_0073;
    const ccomplex_t IT_0075 = e_em*mu_h*m_tau*IT_0000*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0076 = IT_0005*IT_0075;
    const ccomplex_t IT_0077 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0078 = IT_0010*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0081 = IT_0010*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (-0.5)*IT_0074 + 0.5*IT_0076 + (-0.5)*IT_0079 +
       0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0072*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_12 + -IT_0021 + -IT_0022 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = e_em*m_mu*IT_0000*IT_0005;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = IT_0002*IT_0005*IT_0010;
    const ccomplex_t IT_0091 = IT_0073*IT_0090;
    const ccomplex_t IT_0092 = IT_0075*IT_0090;
    const ccomplex_t IT_0093 = IT_0001*IT_0077;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = IT_0001*IT_0080;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (-0.5)*IT_0091 + 0.5*IT_0092 + 0.5*IT_0094 + (
      -0.5)*IT_0096;
    const ccomplex_t IT_0098 = IT_0089*IT_0097;
    const ccomplex_t IT_0099 = pow(m_Z, 2);
    const ccomplex_t IT_0100 = cpow((-2)*s_12 + -IT_0021 + -IT_0022 + IT_0099 
      + -reg_prop, -1);
    const ccomplex_t IT_0101 = IT_0098*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = -IT_0087 + -IT_0102;
    const ccomplex_t IT_0104 = s_13 + s_14;
    const ccomplex_t IT_0105 = m_mu*IT_0104;
    const ccomplex_t IT_0106 = s_13*s_14;
    const ccomplex_t IT_0107 = s_34*IT_0021;
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = IT_0106 + IT_0108;
    const ccomplex_t IT_0110 = IT_0021*IT_0067;
    const ccomplex_t IT_0111 = m_mu*s_14;
    const ccomplex_t IT_0112 = (-2)*IT_0052;
    const ccomplex_t IT_0113 = (-2)*conj(IT_0052);
    const ccomplex_t IT_0114 = (-2)*conj(IT_0066);
    const ccomplex_t IT_0115 = 0.5*IT_0111;
    const ccomplex_t IT_0116 = s_34 + IT_0067;
    const ccomplex_t IT_0117 = 0.5*conj(IT_0065);
    const ccomplex_t IT_0118 = (-2)*IT_0066;
    const ccomplex_t IT_0119 = conj(IT_0103)*IT_0116;
    const ccomplex_t IT_0120 = 2*s_34;
    const ccomplex_t IT_0121 = conj(IT_0038)*IT_0069;
    const ccomplex_t IT_0122 = 2*IT_0038*((conj(IT_0052) + conj(IT_0057))
      *IT_0060 + 2*(conj(IT_0038) + 1./2*conj(IT_0065) + 1./2*conj(IT_0066))
      *IT_0069) + 2*IT_0052*(conj(IT_0038)*IT_0060 + conj(IT_0066)*IT_0070 + 
      -conj(IT_0103)*IT_0105 + 2*conj(IT_0052)*IT_0109) + (conj(IT_0057)*IT_0110
       + conj(IT_0065)*IT_0111)*IT_0112 + 2*IT_0057*(conj(IT_0038)*IT_0060 +
       conj(IT_0065)*IT_0070 + conj(IT_0103)*IT_0105 + 2*conj(IT_0057)*IT_0109 +
       0.5*IT_0110*IT_0113 + IT_0114*IT_0115) + (-2)*IT_0103*((conj(IT_0052) + 
      -conj(IT_0057))*IT_0105 + (-2)*IT_0116*(conj(IT_0103) + 0.25*IT_0114 +
       IT_0117)) + IT_0118*(conj(IT_0065)*IT_0067 + conj(IT_0057)*IT_0111 +
       IT_0119) + 2*IT_0066*(conj(IT_0052)*IT_0070 + 0.5*conj(IT_0066)*IT_0120 +
       IT_0121) + 2*IT_0065*(conj(IT_0057)*IT_0070 + 0.5*IT_0067*IT_0114 +
       IT_0113*IT_0115 + IT_0119 + IT_0117*IT_0120 + IT_0121);
    return create_ccomplex_return(IT_0122);
}

