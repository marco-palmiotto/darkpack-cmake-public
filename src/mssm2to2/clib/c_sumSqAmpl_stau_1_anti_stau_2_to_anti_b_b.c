#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_2_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_2_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
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
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*U_stau_00*conj(U_stau_01) + IT_0001*IT_0002*U_stau_10
      *conj(U_stau_11));
    const ccomplex_t IT_0005 = (-2)*IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0007 = 0.333333333333333*IT_0006;
    const ccomplex_t IT_0008 = pow(m_stau_1, 2);
    const ccomplex_t IT_0009 = pow(m_stau_2, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0008 + -IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = IT_0005*IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = pow(m_b, 2);
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = m_b*s_14;
    const ccomplex_t IT_0018 = -IT_0004;
    const ccomplex_t IT_0019 = m_b*IT_0018;
    const ccomplex_t IT_0020 = m_b*IT_0005;
    const ccomplex_t IT_0021 = IT_0001*IT_0002;
    const ccomplex_t IT_0022 = e_em*IT_0021;
    const ccomplex_t IT_0023 = IT_0000*IT_0003;
    const ccomplex_t IT_0024 = e_em*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0022 + 3*IT_0024);
    const ccomplex_t IT_0026 = (-0.166666666666667)*IT_0025;
    const ccomplex_t IT_0027 = IT_0010*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -(IT_0019 + -IT_0020)*(IT_0012 + -IT_0028);
    const ccomplex_t IT_0030 = IT_0005*IT_0028;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = (-6)*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = m_b*e_em*IT_0003*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = e_em*mu_h*m_tau*IT_0033*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0037 = cos(beta);
    const ccomplex_t IT_0038 = cpow(IT_0037, -1);
    const ccomplex_t IT_0039 = sin(beta);
    const ccomplex_t IT_0040 = IT_0003*IT_0038*IT_0039;
    const ccomplex_t IT_0041 = IT_0036*IT_0040;
    const ccomplex_t IT_0042 = e_em*mu_h*m_tau*IT_0033*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0045 = IT_0037*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0048 = IT_0037*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0041 + 0.5*IT_0043 + 0.5*IT_0046 + (
      -0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0035*IT_0050;
    const ccomplex_t IT_0052 = pow(m_Z, 2);
    const ccomplex_t IT_0053 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + IT_0052 
      + -reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = m_b*e_em*IT_0003*IT_0033*IT_0038*IT_0039;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0003*IT_0036;
    const ccomplex_t IT_0059 = IT_0003*IT_0042;
    const ccomplex_t IT_0060 = IT_0039*IT_0044;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = IT_0039*IT_0047;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (-0.5)*IT_0058 + 0.5*IT_0059 + (-0.5)*IT_0061 +
       0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0057*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = -IT_0055 + -IT_0068;
    const ccomplex_t IT_0070 = s_13 + s_14;
    const ccomplex_t IT_0071 = m_b*IT_0070;
    const ccomplex_t IT_0072 = s_34 + IT_0015;
    const ccomplex_t IT_0073 = -IT_0029;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = (-6)*conj(IT_0075);
    const ccomplex_t IT_0077 = 12*s_34;
    const ccomplex_t IT_0078 = cos(alpha);
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*m_b*e_em*IT_0003*IT_0033
      *IT_0038*IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = mu_h*m_tau*IT_0033;
    const ccomplex_t IT_0082 = sin(alpha);
    const ccomplex_t IT_0083 = -IT_0039*IT_0082;
    const ccomplex_t IT_0084 = IT_0037*IT_0078;
    const ccomplex_t IT_0085 = IT_0083 + IT_0084;
    const ccomplex_t IT_0086 = cpow(IT_0000, -2);
    const ccomplex_t IT_0087 = pow(m_tau, 2);
    const ccomplex_t IT_0088 = (-2)*U_stau_00*(((-0.5)*m_W*IT_0085*(IT_0003 + 
      -IT_0002*IT_0086) + IT_0003*IT_0033*IT_0038*IT_0078*IT_0087)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0038*IT_0081*IT_0082*conj(U_stau_11)) +
       U_stau_10*(IT_0003*IT_0038*IT_0081*IT_0082*conj(U_stau_01) + 2*(m_W
      *IT_0002*IT_0085*IT_0086 + -IT_0003*IT_0033*IT_0038*IT_0078*IT_0087)*conj
      (U_stau_11));
    const ccomplex_t IT_0089 = (0 + _Complex_I*-1)*(e_em*IT_0088 + (
      -1.4142135623731)*A_tau*IT_0078*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = IT_0080*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*m_b*e_em*IT_0003*IT_0033
      *IT_0038*IT_0082;
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0037*IT_0082;
    const ccomplex_t IT_0098 = IT_0039*IT_0078;
    const ccomplex_t IT_0099 = IT_0097 + IT_0098;
    const ccomplex_t IT_0100 = 2*U_stau_00*((IT_0003*IT_0033*IT_0038*IT_0082
      *IT_0087 + (-0.5)*m_W*(IT_0003 + -IT_0002*IT_0086)*IT_0099)*conj(U_stau_01
      ) + 1./2*IT_0003*IT_0038*IT_0078*IT_0081*conj(U_stau_11)) + U_stau_10*
      (IT_0003*IT_0038*IT_0078*IT_0081*conj(U_stau_01) + 2*(IT_0003*IT_0033
      *IT_0038*IT_0082*IT_0087 + -m_W*IT_0002*IT_0086*IT_0099)*conj(U_stau_11));
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(e_em*IT_0100 +
       1.4142135623731*A_tau*IT_0082*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = IT_0096*IT_0102;
    const ccomplex_t IT_0104 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0008 + -IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = -IT_0094 + -IT_0106;
    const ccomplex_t IT_0108 = -s_14;
    const ccomplex_t IT_0109 = s_13 + IT_0108;
    const ccomplex_t IT_0110 = m_b*IT_0109;
    const ccomplex_t IT_0111 = -IT_0015;
    const ccomplex_t IT_0112 = s_34 + IT_0111;
    const ccomplex_t IT_0113 = m_b*s_13;
    const ccomplex_t IT_0114 = s_34*IT_0008;
    const ccomplex_t IT_0115 = (-6)*IT_0114;
    const ccomplex_t IT_0116 = s_13*s_14;
    const ccomplex_t IT_0117 = 12*IT_0116;
    const ccomplex_t IT_0118 = IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = conj(IT_0107)*IT_0110;
    const ccomplex_t IT_0120 = (-6)*conj(IT_0031);
    const ccomplex_t IT_0121 = 0.166666666666667*IT_0017;
    const ccomplex_t IT_0122 = 6*s_34;
    const ccomplex_t IT_0123 = conj(IT_0107)*IT_0112;
    const ccomplex_t IT_0124 = (-6)*IT_0075;
    const ccomplex_t IT_0125 = conj(IT_0069)*IT_0072;
    const ccomplex_t IT_0126 = (conj(IT_0014)*IT_0016 + IT_0017*conj(IT_0029))
      *IT_0032 + 6*IT_0069*((conj(IT_0014) + -conj(IT_0031))*IT_0071 + IT_0072*
      (conj(IT_0029) + 0.166666666666667*IT_0076) + 0.166666666666667*conj
      (IT_0069)*(12*IT_0015 + IT_0077)) + 6*IT_0107*(0.166666666666667*((-12)
      *IT_0015 + IT_0077)*conj(IT_0107) + (conj(IT_0014) + conj(IT_0031))
      *IT_0110 + (conj(IT_0029) + conj(IT_0075))*IT_0112) + (-6)*IT_0031*(conj
      (IT_0069)*IT_0071 + -conj(IT_0075)*IT_0113 + (-0.166666666666667)*conj
      (IT_0031)*IT_0118 + -IT_0119) + 6*IT_0014*(conj(IT_0069)*IT_0071 + conj
      (IT_0029)*IT_0113 + 0.166666666666667*conj(IT_0014)*IT_0118 + IT_0119 +
       0.166666666666667*IT_0016*IT_0120 + IT_0076*IT_0121) + 6*IT_0075*(conj
      (IT_0031)*IT_0113 + 0.166666666666667*conj(IT_0075)*IT_0122 + IT_0123) +
       IT_0124*(conj(IT_0014)*IT_0017 + IT_0015*conj(IT_0029) + IT_0125) + 6
      *IT_0029*(0.166666666666667*IT_0015*IT_0076 + conj(IT_0014)*IT_0113 +
       IT_0120*IT_0121 + 0.166666666666667*conj(IT_0029)*IT_0122 + IT_0123 +
       IT_0125);
    return create_ccomplex_return(IT_0126);
}

