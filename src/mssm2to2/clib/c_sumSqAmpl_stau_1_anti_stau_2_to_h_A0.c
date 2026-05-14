#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_2_to_h_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_2_to_h_A0(
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
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
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
    const ccomplex_t IT_0000 = pow(m_stau_1, 2);
    const ccomplex_t IT_0001 = pow(m_stau_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0000 + -IT_0001 + -reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cos(alpha);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = e_em*IT_0004*(IT_0005*IT_0006 + IT_0007*IT_0008
      )*IT_0010;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0004*IT_0009 +
       -IT_0003*IT_0010)*U_stau_00*conj(U_stau_01) + IT_0004*IT_0009*U_stau_10
      *conj(U_stau_11));
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = s_13 + s_23;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (-0.5)*IT_0011;
    const ccomplex_t IT_0019 = (-2)*IT_0013;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = s_12 + IT_0001;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = IT_0017 + IT_0022;
    const ccomplex_t IT_0024 = IT_0012*IT_0019;
    const ccomplex_t IT_0025 = s_23*IT_0024;
    const ccomplex_t IT_0026 = IT_0014*IT_0018;
    const ccomplex_t IT_0027 = 0.5*IT_0000 + 0.5*IT_0001;
    const ccomplex_t IT_0028 = s_12 + IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = -IT_0025 + (-2)*IT_0029;
    const ccomplex_t IT_0031 = IT_0023 + IT_0030;
    const ccomplex_t IT_0032 = IT_0002*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_12 + -IT_0000 + -IT_0001 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0035 = pow(m_W, -1);
    const ccomplex_t IT_0036 = e_em*mu_h*m_tau*IT_0035*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0037 = IT_0010*IT_0036;
    const ccomplex_t IT_0038 = e_em*mu_h*m_tau*IT_0035*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0039 = IT_0010*IT_0038;
    const ccomplex_t IT_0040 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0041 = IT_0007*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0044 = IT_0007*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0037 + 0.5*IT_0039 + (-0.5)*IT_0042 +
       0.5*IT_0045;
    const ccomplex_t IT_0047 = cpow(IT_0003, -2);
    const ccomplex_t IT_0048 = IT_0009*IT_0047;
    const ccomplex_t IT_0049 = IT_0010 + IT_0048;
    const ccomplex_t IT_0050 = cpow(IT_0007, 3);
    const ccomplex_t IT_0051 = cpow(IT_0007, 2);
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*m_W*e_em*IT_0049*(IT_0006
      *IT_0050 + -IT_0005*(IT_0005*(IT_0006*IT_0007 + IT_0005*IT_0008) + 
      -IT_0008*IT_0051));
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0046*IT_0053;
    const ccomplex_t IT_0055 = IT_0034*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = pow(m_h, 2);
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0057 + reg_prop + (0 +
       _Complex_I*1)*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0059 = mu_h*m_tau*IT_0035;
    const ccomplex_t IT_0060 = cpow(IT_0005, -1);
    const ccomplex_t IT_0061 = IT_0005*IT_0008;
    const ccomplex_t IT_0062 = IT_0006*IT_0007;
    const ccomplex_t IT_0063 = IT_0061 + IT_0062;
    const ccomplex_t IT_0064 = pow(m_tau, 2);
    const ccomplex_t IT_0065 = 2*U_stau_00*(((-0.5)*m_W*(IT_0010 + -IT_0009
      *IT_0047)*IT_0063 + IT_0008*IT_0010*IT_0035*IT_0060*IT_0064)*conj
      (U_stau_00) + 1./2*IT_0006*IT_0010*IT_0059*IT_0060*conj(U_stau_10)) +
       U_stau_10*(IT_0006*IT_0010*IT_0059*IT_0060*conj(U_stau_00) + (-2)*(m_W
      *IT_0009*IT_0047*IT_0063 + -IT_0008*IT_0010*IT_0035*IT_0060*IT_0064)*conj
      (U_stau_10));
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(e_em*IT_0065 +
       1.4142135623731*A_tau*IT_0008*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0046*IT_0067;
    const ccomplex_t IT_0069 = IT_0058*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_23 + IT_0057 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0072 = 2*U_stau_01*(((-0.5)*m_W*(IT_0010 + -IT_0009
      *IT_0047)*IT_0063 + IT_0008*IT_0010*IT_0035*IT_0060*IT_0064)*conj
      (U_stau_01) + 1./2*IT_0006*IT_0010*IT_0059*IT_0060*conj(U_stau_11)) +
       U_stau_11*(IT_0006*IT_0010*IT_0059*IT_0060*conj(U_stau_01) + (-2)*(m_W
      *IT_0009*IT_0047*IT_0063 + -IT_0008*IT_0010*IT_0035*IT_0060*IT_0064)*conj
      (U_stau_11));
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(e_em*IT_0072 +
       1.4142135623731*A_tau*IT_0008*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0046*IT_0074;
    const ccomplex_t IT_0076 = IT_0071*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0001 + IT_0057 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0079 = e_em*mu_h*m_tau*IT_0035*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0080 = IT_0010*IT_0079;
    const ccomplex_t IT_0081 = e_em*mu_h*m_tau*IT_0035*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0082 = IT_0010*IT_0081;
    const ccomplex_t IT_0083 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0084 = IT_0007*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0087 = IT_0007*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (-0.5)*IT_0080 + 0.5*IT_0082 + (-0.5)*IT_0085 +
       0.5*IT_0088;
    const ccomplex_t IT_0090 = 2*U_stau_00*(((-0.5)*m_W*(IT_0010 + -IT_0009
      *IT_0047)*IT_0063 + IT_0008*IT_0010*IT_0035*IT_0060*IT_0064)*conj
      (U_stau_01) + 1./2*IT_0006*IT_0010*IT_0059*IT_0060*conj(U_stau_11)) +
       U_stau_10*(IT_0006*IT_0010*IT_0059*IT_0060*conj(U_stau_01) + (-2)*(m_W
      *IT_0009*IT_0047*IT_0063 + -IT_0008*IT_0010*IT_0035*IT_0060*IT_0064)*conj
      (U_stau_11));
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(e_em*IT_0090 +
       1.4142135623731*A_tau*IT_0008*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0089*IT_0092;
    const ccomplex_t IT_0094 = IT_0078*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_13 + IT_0000 + IT_0057 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0097 = e_em*mu_h*m_tau*IT_0035*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0098 = IT_0010*IT_0097;
    const ccomplex_t IT_0099 = e_em*mu_h*m_tau*IT_0035*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0100 = IT_0010*IT_0099;
    const ccomplex_t IT_0101 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0102 = IT_0007*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0105 = IT_0007*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (-0.5)*IT_0098 + 0.5*IT_0100 + (-0.5)*IT_0103 +
       0.5*IT_0106;
    const ccomplex_t IT_0108 = IT_0092*IT_0107;
    const ccomplex_t IT_0109 = IT_0096*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*m_W*e_em*IT_0005*IT_0007*
      (IT_0006*IT_0007 + IT_0005*IT_0008)*IT_0049;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0007*IT_0010*IT_0060;
    const ccomplex_t IT_0114 = IT_0036*IT_0113;
    const ccomplex_t IT_0115 = IT_0038*IT_0113;
    const ccomplex_t IT_0116 = IT_0005*IT_0040;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = IT_0005*IT_0043;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = (-0.5)*IT_0114 + 0.5*IT_0115 + 0.5*IT_0117 + (
      -0.5)*IT_0119;
    const ccomplex_t IT_0121 = IT_0112*IT_0120;
    const ccomplex_t IT_0122 = pow(m_Z, 2);
    const ccomplex_t IT_0123 = cpow((-2)*s_12 + -IT_0000 + -IT_0001 + IT_0122 
      + -reg_prop, -1);
    const ccomplex_t IT_0124 = IT_0121*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = IT_0033 + -IT_0056 + IT_0070 + IT_0077 +
       IT_0095 + IT_0110 + -IT_0125;
    return create_ccomplex_return(IT_0126*conj(IT_0126));
}

