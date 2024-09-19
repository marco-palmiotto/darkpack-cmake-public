#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_anti_tau_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
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
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_00*conj(U_sb_01) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_11));
    const ccomplex_t IT_0005 = (-0.666666666666667)*IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0007 = pow(m_sb_1, 2);
    const ccomplex_t IT_0008 = pow(m_sb_2, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0007 + -IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = pow(m_tau, 2);
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = s_14*m_tau;
    const ccomplex_t IT_0017 = (-0.333333333333333)*IT_0004;
    const ccomplex_t IT_0018 = m_tau*IT_0017;
    const ccomplex_t IT_0019 = m_tau*IT_0005;
    const ccomplex_t IT_0020 = IT_0001*IT_0002;
    const ccomplex_t IT_0021 = e_em*IT_0020;
    const ccomplex_t IT_0022 = IT_0000*IT_0003;
    const ccomplex_t IT_0023 = e_em*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0021 + -IT_0023);
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (IT_0018 + -IT_0019)*(IT_0011 + -IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0005*IT_0027;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = (-6)*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = cos(beta);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = sin(beta);
    const ccomplex_t IT_0037 = e_em*m_tau*IT_0003*IT_0033*IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = m_b*e_em*mu_h*IT_0033*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0040 = IT_0003*IT_0039;
    const ccomplex_t IT_0041 = m_b*e_em*mu_h*IT_0033*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0042 = IT_0003*IT_0041;
    const ccomplex_t IT_0043 = cpow(V_tb, 2);
    const ccomplex_t IT_0044 = A_b*IT_0043*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0045 = IT_0036*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = A_b*IT_0043*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0048 = IT_0036*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0040 + 0.5*IT_0042 + (-0.5)*IT_0046 +
       0.5*IT_0049;
    const ccomplex_t IT_0051 = IT_0038*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = e_em*m_tau*IT_0003*IT_0033;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = IT_0003*IT_0035*IT_0036;
    const ccomplex_t IT_0058 = IT_0039*IT_0057;
    const ccomplex_t IT_0059 = IT_0041*IT_0057;
    const ccomplex_t IT_0060 = IT_0034*IT_0044;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = IT_0034*IT_0047;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (-0.5)*IT_0058 + 0.5*IT_0059 + 0.5*IT_0061 + (
      -0.5)*IT_0063;
    const ccomplex_t IT_0065 = IT_0056*IT_0064;
    const ccomplex_t IT_0066 = pow(m_Z, 2);
    const ccomplex_t IT_0067 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + IT_0066 
      + -reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = -IT_0054 + -IT_0069;
    const ccomplex_t IT_0071 = s_13 + s_14;
    const ccomplex_t IT_0072 = m_tau*IT_0071;
    const ccomplex_t IT_0073 = s_34 + IT_0014;
    const ccomplex_t IT_0074 = -IT_0029;
    const ccomplex_t IT_0075 = (-6)*conj(IT_0074);
    const ccomplex_t IT_0076 = 12*s_34;
    const ccomplex_t IT_0077 = cos(alpha);
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*e_em*m_tau*IT_0003*IT_0033
      *IT_0035*IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = m_b*mu_h*IT_0033;
    const ccomplex_t IT_0081 = sin(alpha);
    const ccomplex_t IT_0082 = -IT_0036*IT_0081;
    const ccomplex_t IT_0083 = IT_0034*IT_0077;
    const ccomplex_t IT_0084 = IT_0082 + IT_0083;
    const ccomplex_t IT_0085 = cpow(IT_0000, -2);
    const ccomplex_t IT_0086 = pow(m_b, 2);
    const ccomplex_t IT_0087 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0084*(IT_0003 + 1.
      /3*IT_0002*IT_0085) + IT_0003*IT_0033*IT_0035*IT_0077*IT_0086)*conj
      (U_sb_01) + -1./2*IT_0003*IT_0035*IT_0080*IT_0081*conj(U_sb_11)) + U_sb_10
      *(IT_0003*IT_0035*IT_0080*IT_0081*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0002*IT_0084*IT_0085 + (-3)*IT_0003*IT_0033*IT_0035*IT_0077*IT_0086)
      *conj(U_sb_11));
    const ccomplex_t IT_0088 = (0 + _Complex_I*-1)*(e_em*IT_0087 + (
      -1.4142135623731)*A_b*IT_0043*IT_0077*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0079*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*e_em*m_tau*IT_0003*IT_0033
      *IT_0035*IT_0081;
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = IT_0034*IT_0081;
    const ccomplex_t IT_0097 = IT_0036*IT_0077;
    const ccomplex_t IT_0098 = IT_0096 + IT_0097;
    const ccomplex_t IT_0099 = 2*U_sb_00*((IT_0003*IT_0033*IT_0035*IT_0081
      *IT_0086 + (-0.5)*m_W*(IT_0003 + 1./3*IT_0002*IT_0085)*IT_0098)*conj
      (U_sb_01) + 1./2*IT_0003*IT_0035*IT_0077*IT_0080*conj(U_sb_11)) + U_sb_10*
      (IT_0003*IT_0035*IT_0077*IT_0080*conj(U_sb_01) + 2*(IT_0003*IT_0033
      *IT_0035*IT_0081*IT_0086 + -1./3*m_W*IT_0002*IT_0085*IT_0098)*conj(U_sb_11
      ));
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(e_em*IT_0099 +
       1.4142135623731*A_b*IT_0043*IT_0081*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = IT_0095*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0007 + -IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = -IT_0093 + -IT_0105;
    const ccomplex_t IT_0107 = -s_14;
    const ccomplex_t IT_0108 = s_13 + IT_0107;
    const ccomplex_t IT_0109 = m_tau*IT_0108;
    const ccomplex_t IT_0110 = -IT_0014;
    const ccomplex_t IT_0111 = s_34 + IT_0110;
    const ccomplex_t IT_0112 = (-6)*IT_0074;
    const ccomplex_t IT_0113 = conj(IT_0070)*IT_0073;
    const ccomplex_t IT_0114 = s_13*m_tau;
    const ccomplex_t IT_0115 = s_34*IT_0007;
    const ccomplex_t IT_0116 = (-6)*IT_0115;
    const ccomplex_t IT_0117 = s_13*s_14;
    const ccomplex_t IT_0118 = 12*IT_0117;
    const ccomplex_t IT_0119 = IT_0116 + IT_0118;
    const ccomplex_t IT_0120 = conj(IT_0106)*IT_0109;
    const ccomplex_t IT_0121 = (-6)*conj(IT_0031);
    const ccomplex_t IT_0122 = 0.166666666666667*IT_0016;
    const ccomplex_t IT_0123 = 6*s_34;
    const ccomplex_t IT_0124 = conj(IT_0106)*IT_0111;
    const ccomplex_t IT_0125 = (conj(IT_0013)*IT_0015 + IT_0016*conj(IT_0029))
      *IT_0032 + 6*IT_0070*((conj(IT_0013) + -conj(IT_0031))*IT_0072 + IT_0073*
      (conj(IT_0029) + 0.166666666666667*IT_0075) + 0.166666666666667*conj
      (IT_0070)*(12*IT_0014 + IT_0076)) + 6*IT_0106*(0.166666666666667*((-12)
      *IT_0014 + IT_0076)*conj(IT_0106) + (conj(IT_0013) + conj(IT_0031))
      *IT_0109 + (conj(IT_0029) + conj(IT_0074))*IT_0111) + IT_0112*(conj
      (IT_0013)*IT_0016 + IT_0014*conj(IT_0029) + IT_0113) + (-6)*IT_0031*(conj
      (IT_0070)*IT_0072 + -conj(IT_0074)*IT_0114 + (-0.166666666666667)*conj
      (IT_0031)*IT_0119 + -IT_0120) + 6*IT_0013*(conj(IT_0070)*IT_0072 + conj
      (IT_0029)*IT_0114 + 0.166666666666667*conj(IT_0013)*IT_0119 + IT_0120 +
       0.166666666666667*IT_0015*IT_0121 + IT_0075*IT_0122) + 6*IT_0029*
      (0.166666666666667*IT_0014*IT_0075 + IT_0113 + conj(IT_0013)*IT_0114 +
       IT_0121*IT_0122 + 0.166666666666667*conj(IT_0029)*IT_0123 + IT_0124) + 6
      *IT_0074*(conj(IT_0031)*IT_0114 + 0.166666666666667*conj(IT_0074)*IT_0123 
      + IT_0124);
    return create_ccomplex_return(IT_0125);
}

