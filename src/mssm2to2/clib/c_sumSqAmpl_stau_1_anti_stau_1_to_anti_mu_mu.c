#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_1_to_anti_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_1_to_anti_mu_mu(
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
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*U_stau_00*conj(U_stau_00) + IT_0001*IT_0002*U_stau_10
      *conj(U_stau_10));
    const ccomplex_t IT_0005 = (-2)*IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0002;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0003;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0007 + -IT_0009);
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = pow(m_stau_1, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0005*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = -IT_0004;
    const ccomplex_t IT_0019 = m_mu*IT_0018;
    const ccomplex_t IT_0020 = m_mu*IT_0005;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0022 = IT_0013*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = (IT_0019 + -IT_0020)*(IT_0015 + -IT_0023);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = m_mu*s_13;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = e_em*m_mu*IT_0003*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = e_em*mu_h*m_tau*IT_0029*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0033 = cos(beta);
    const ccomplex_t IT_0034 = cpow(IT_0033, -1);
    const ccomplex_t IT_0035 = sin(beta);
    const ccomplex_t IT_0036 = IT_0003*IT_0034*IT_0035;
    const ccomplex_t IT_0037 = IT_0032*IT_0036;
    const ccomplex_t IT_0038 = e_em*mu_h*m_tau*IT_0029*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0041 = IT_0033*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0044 = IT_0033*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0037 + 0.5*IT_0039 + 0.5*IT_0042 + (
      -0.5)*IT_0045;
    const ccomplex_t IT_0047 = IT_0031*IT_0046;
    const ccomplex_t IT_0048 = pow(m_Z, 2);
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + (-2)*IT_0012 + IT_0048 + 
      -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = e_em*m_mu*IT_0003*IT_0029*IT_0034*IT_0035;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = IT_0003*IT_0032;
    const ccomplex_t IT_0055 = IT_0003*IT_0038;
    const ccomplex_t IT_0056 = IT_0035*IT_0040;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = IT_0035*IT_0043;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (-0.5)*IT_0054 + 0.5*IT_0055 + (-0.5)*IT_0057 +
       0.5*IT_0059;
    const ccomplex_t IT_0061 = IT_0053*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_12 + (-2)*IT_0012 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = -IT_0051 + -IT_0064;
    const ccomplex_t IT_0066 = s_13 + s_14;
    const ccomplex_t IT_0067 = m_mu*IT_0066;
    const ccomplex_t IT_0068 = s_13*s_14;
    const ccomplex_t IT_0069 = s_34*IT_0012;
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = IT_0068 + IT_0070;
    const ccomplex_t IT_0072 = 2*s_34;
    const ccomplex_t IT_0073 = pow(m_mu, 2);
    const ccomplex_t IT_0074 = IT_0012*IT_0073;
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = IT_0068 + IT_0070 + IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = 2*IT_0077;
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = cpow(s_12 + IT_0012 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = cpow(IT_0084, 2);
    const ccomplex_t IT_0086 = IT_0005*IT_0023;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = conj(IT_0017) + conj(IT_0087);
    const ccomplex_t IT_0089 = 4*IT_0084;
    const ccomplex_t IT_0090 = (-4)*IT_0084;
    const ccomplex_t IT_0091 = m_mu*s_14;
    const ccomplex_t IT_0092 = (-2)*IT_0017;
    const ccomplex_t IT_0093 = -s_14;
    const ccomplex_t IT_0094 = s_13 + IT_0093;
    const ccomplex_t IT_0095 = m_mu*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0097 = sin(alpha);
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*e_em*m_mu*IT_0003*IT_0029
      *IT_0034*IT_0097;
    const ccomplex_t IT_0099 = 0.5*IT_0098;
    const ccomplex_t IT_0100 = mu_h*m_tau*IT_0029;
    const ccomplex_t IT_0101 = cos(alpha);
    const ccomplex_t IT_0102 = IT_0033*IT_0097;
    const ccomplex_t IT_0103 = IT_0035*IT_0101;
    const ccomplex_t IT_0104 = IT_0102 + IT_0103;
    const ccomplex_t IT_0105 = cpow(IT_0000, -2);
    const ccomplex_t IT_0106 = pow(m_tau, 2);
    const ccomplex_t IT_0107 = 2*U_stau_00*(((-0.5)*m_W*IT_0104*(IT_0003 + 
      -IT_0002*IT_0105) + IT_0003*IT_0029*IT_0034*IT_0097*IT_0106)*conj
      (U_stau_00) + 1./2*IT_0003*IT_0034*IT_0100*IT_0101*conj(U_stau_10)) +
       U_stau_10*(IT_0003*IT_0034*IT_0100*IT_0101*conj(U_stau_00) + (-2)*(m_W
      *IT_0002*IT_0104*IT_0105 + -IT_0003*IT_0029*IT_0034*IT_0097*IT_0106)*conj
      (U_stau_10));
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(e_em*IT_0107 +
       1.4142135623731*A_tau*IT_0097*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = IT_0099*IT_0109;
    const ccomplex_t IT_0111 = IT_0096*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = cpow((-2)*s_12 + (-2)*IT_0012 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*e_em*m_mu*IT_0003*IT_0029
      *IT_0034*IT_0101;
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = -IT_0035*IT_0097;
    const ccomplex_t IT_0117 = IT_0033*IT_0101;
    const ccomplex_t IT_0118 = IT_0116 + IT_0117;
    const ccomplex_t IT_0119 = (-2)*U_stau_00*((IT_0003*IT_0029*IT_0034
      *IT_0101*IT_0106 + (-0.5)*m_W*(IT_0003 + -IT_0002*IT_0105)*IT_0118)*conj
      (U_stau_00) + -1./2*IT_0003*IT_0034*IT_0097*IT_0100*conj(U_stau_10)) +
       U_stau_10*(IT_0003*IT_0034*IT_0097*IT_0100*conj(U_stau_00) + (-2)*
      (IT_0003*IT_0029*IT_0034*IT_0101*IT_0106 + -m_W*IT_0002*IT_0105*IT_0118)
      *conj(U_stau_10));
    const ccomplex_t IT_0120 = (0 + _Complex_I*-1)*(e_em*IT_0119 + (
      -1.4142135623731)*A_tau*IT_0101*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = IT_0115*IT_0121;
    const ccomplex_t IT_0123 = IT_0113*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = -IT_0112 + -IT_0124;
    const ccomplex_t IT_0126 = (-2)*IT_0027;
    const ccomplex_t IT_0127 = (-2)*conj(IT_0017);
    const ccomplex_t IT_0128 = (-2)*conj(IT_0027);
    const ccomplex_t IT_0129 = 0.5*IT_0091;
    const ccomplex_t IT_0130 = -IT_0073;
    const ccomplex_t IT_0131 = s_34 + IT_0130;
    const ccomplex_t IT_0132 = 0.5*conj(IT_0026);
    const ccomplex_t IT_0133 = s_34 + IT_0073;
    const ccomplex_t IT_0134 = conj(IT_0065)*IT_0133;
    const ccomplex_t IT_0135 = 2*IT_0017*(conj(IT_0027)*IT_0028 + -conj
      (IT_0065)*IT_0067 + 2*conj(IT_0017)*IT_0071) + 2*IT_0027*(conj(IT_0017)
      *IT_0028 + 0.5*conj(IT_0027)*IT_0072) + (-8)*IT_0076*(IT_0085 + (-0.125)
      *IT_0088*IT_0089 + (-0.125)*IT_0017*IT_0090) + (IT_0074*conj(IT_0087) +
       conj(IT_0026)*IT_0091)*IT_0092 + 2*IT_0095*(IT_0088*IT_0125 + (IT_0017 +
       IT_0087)*conj(IT_0125) + (-2)*IT_0084*(0.5*IT_0026 + (-0.5)*conj(IT_0026)
       + (-0.5)*conj(IT_0027) + IT_0125 + -conj(IT_0125) + (-0.25)*IT_0126)) + 2
      *IT_0087*(conj(IT_0026)*IT_0028 + conj(IT_0065)*IT_0067 + 2*IT_0071*conj
      (IT_0087) + 0.5*IT_0076*IT_0090 + 0.5*IT_0074*IT_0127 + IT_0128*IT_0129) +
       2*IT_0131*(IT_0027*conj(IT_0125) + 2*IT_0125*(0.5*conj(IT_0027) + conj
      (IT_0125) + IT_0132)) + (-2)*IT_0065*(IT_0067*(conj(IT_0017) + -conj
      (IT_0087)) + (-2)*(conj(IT_0065) + 0.25*IT_0128 + IT_0132)*IT_0133) +
       IT_0126*(conj(IT_0026)*IT_0073 + conj(IT_0087)*IT_0091 + IT_0134) + 2
      *IT_0026*(IT_0028*conj(IT_0087) + 0.5*IT_0073*IT_0128 + IT_0127*IT_0129 +
       conj(IT_0125)*IT_0131 + IT_0072*IT_0132 + IT_0134);
    return create_ccomplex_return(IT_0135);
}

