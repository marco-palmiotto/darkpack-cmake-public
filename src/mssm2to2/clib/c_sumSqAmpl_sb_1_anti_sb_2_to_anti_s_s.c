#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_anti_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
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
    const creal_t alpha = param->alpha;
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
    const ccomplex_t IT_0006 = IT_0001*IT_0002;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0003;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0007 + 3*IT_0009);
    const ccomplex_t IT_0011 = (-0.166666666666667)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_sb_1, 2);
    const ccomplex_t IT_0013 = pow(m_sb_2, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0012 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0005*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = (-0.333333333333333)*IT_0004;
    const ccomplex_t IT_0020 = m_s*IT_0019;
    const ccomplex_t IT_0021 = m_s*IT_0005;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0023 = 0.333333333333333*IT_0022;
    const ccomplex_t IT_0024 = IT_0014*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (IT_0020 + -IT_0021)*(IT_0016 + -IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = m_s*s_13;
    const ccomplex_t IT_0031 = pow(m_W, -1);
    const ccomplex_t IT_0032 = cos(beta);
    const ccomplex_t IT_0033 = cpow(IT_0032, -1);
    const ccomplex_t IT_0034 = sin(beta);
    const ccomplex_t IT_0035 = m_s*e_em*IT_0003*IT_0031*IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = m_b*e_em*mu_h*IT_0031*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = m_b*e_em*mu_h*IT_0031*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0040 = IT_0003*IT_0039;
    const ccomplex_t IT_0041 = cpow(V_tb, 2);
    const ccomplex_t IT_0042 = A_b*IT_0041*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0043 = IT_0034*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = A_b*IT_0041*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0046 = IT_0034*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (-0.5)*IT_0038 + 0.5*IT_0040 + (-0.5)*IT_0044 +
       0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0036*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = m_s*e_em*IT_0003*IT_0031;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0003*IT_0033*IT_0034;
    const ccomplex_t IT_0056 = IT_0037*IT_0055;
    const ccomplex_t IT_0057 = IT_0039*IT_0055;
    const ccomplex_t IT_0058 = IT_0032*IT_0042;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = IT_0032*IT_0045;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (-0.5)*IT_0056 + 0.5*IT_0057 + 0.5*IT_0059 + (
      -0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0054*IT_0062;
    const ccomplex_t IT_0064 = pow(m_Z, 2);
    const ccomplex_t IT_0065 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + IT_0064 
      + -reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0052 + -IT_0067;
    const ccomplex_t IT_0069 = s_13 + s_14;
    const ccomplex_t IT_0070 = m_s*IT_0069;
    const ccomplex_t IT_0071 = sin(alpha);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_s*e_em*IT_0003*IT_0031
      *IT_0033*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = m_b*mu_h*IT_0031;
    const ccomplex_t IT_0075 = cos(alpha);
    const ccomplex_t IT_0076 = IT_0032*IT_0071;
    const ccomplex_t IT_0077 = IT_0034*IT_0075;
    const ccomplex_t IT_0078 = IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = cpow(IT_0000, -2);
    const ccomplex_t IT_0080 = pow(m_b, 2);
    const ccomplex_t IT_0081 = 2*U_sb_00*(((-0.5)*m_W*IT_0078*(IT_0003 + 1./3
      *IT_0002*IT_0079) + IT_0003*IT_0031*IT_0033*IT_0071*IT_0080)*conj(U_sb_01)
       + 1./2*IT_0003*IT_0033*IT_0074*IT_0075*conj(U_sb_11)) + U_sb_10*(IT_0003
      *IT_0033*IT_0074*IT_0075*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0002
      *IT_0078*IT_0079 + (-3)*IT_0003*IT_0031*IT_0033*IT_0071*IT_0080)*conj
      (U_sb_11));
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(e_em*IT_0081 +
       1.4142135623731*A_b*IT_0041*IT_0071*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0073*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0012 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*m_s*e_em*IT_0003*IT_0031
      *IT_0033*IT_0075;
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = -IT_0034*IT_0071;
    const ccomplex_t IT_0091 = IT_0032*IT_0075;
    const ccomplex_t IT_0092 = IT_0090 + IT_0091;
    const ccomplex_t IT_0093 = (-2)*U_sb_00*((IT_0003*IT_0031*IT_0033*IT_0075
      *IT_0080 + (-0.5)*m_W*(IT_0003 + 1./3*IT_0002*IT_0079)*IT_0092)*conj
      (U_sb_01) + -1./2*IT_0003*IT_0033*IT_0071*IT_0074*conj(U_sb_11)) + U_sb_10
      *(IT_0003*IT_0033*IT_0071*IT_0074*conj(U_sb_01) + (-2)*(IT_0003*IT_0031
      *IT_0033*IT_0075*IT_0080 + -1./3*m_W*IT_0002*IT_0079*IT_0092)*conj(U_sb_11
      ));
    const ccomplex_t IT_0094 = (0 + _Complex_I*-1)*(e_em*IT_0093 + (
      -1.4142135623731)*A_b*IT_0041*IT_0075*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0089*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = -IT_0087 + -IT_0099;
    const ccomplex_t IT_0101 = -s_14;
    const ccomplex_t IT_0102 = s_13 + IT_0101;
    const ccomplex_t IT_0103 = m_s*IT_0102;
    const ccomplex_t IT_0104 = s_34*IT_0012;
    const ccomplex_t IT_0105 = (-18)*IT_0104;
    const ccomplex_t IT_0106 = s_13*s_14;
    const ccomplex_t IT_0107 = 36*IT_0106;
    const ccomplex_t IT_0108 = IT_0105 + IT_0107;
    const ccomplex_t IT_0109 = IT_0005*IT_0025;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = m_s*s_14;
    const ccomplex_t IT_0112 = (-18)*IT_0029;
    const ccomplex_t IT_0113 = pow(m_s, 2);
    const ccomplex_t IT_0114 = IT_0012*IT_0113;
    const ccomplex_t IT_0115 = (-18)*IT_0018;
    const ccomplex_t IT_0116 = (-18)*conj(IT_0029);
    const ccomplex_t IT_0117 = (-18)*conj(IT_0018);
    const ccomplex_t IT_0118 = 18*s_34;
    const ccomplex_t IT_0119 = 18*IT_0113 + IT_0118;
    const ccomplex_t IT_0120 = (-18)*IT_0113;
    const ccomplex_t IT_0121 = (-18)*s_34 + IT_0120;
    const ccomplex_t IT_0122 = 18*IT_0110;
    const ccomplex_t IT_0123 = IT_0118 + IT_0120;
    const ccomplex_t IT_0124 = 18*IT_0018;
    const ccomplex_t IT_0125 = IT_0100*IT_0123;
    const ccomplex_t IT_0126 = conj(IT_0018)*(18*IT_0029*IT_0030 + (-18)
      *IT_0068*IT_0070 + 18*IT_0100*IT_0103 + IT_0018*IT_0108) + conj(IT_0110)*
      (18*IT_0026*IT_0030 + 18*IT_0068*IT_0070 + 18*IT_0100*IT_0103 + IT_0108
      *IT_0110 + IT_0111*IT_0112 + IT_0114*IT_0115) + IT_0026*(IT_0113*IT_0116 +
       IT_0111*IT_0117) + IT_0110*(IT_0111*IT_0116 + IT_0114*IT_0117) + conj
      (IT_0068)*(IT_0068*(36*s_34 + 36*IT_0113) + IT_0026*IT_0119 + IT_0029
      *IT_0121 + IT_0070*((-18)*IT_0018 + IT_0122)) + conj(IT_0100)*(IT_0100*(36
      *s_34 + (-36)*IT_0113) + (IT_0026 + IT_0029)*IT_0123 + IT_0103*(IT_0122 +
       IT_0124)) + conj(IT_0026)*(IT_0112*IT_0113 + IT_0111*IT_0115 + IT_0026
      *IT_0118 + IT_0068*IT_0119 + IT_0030*IT_0122 + IT_0125) + conj(IT_0029)*
      (IT_0029*IT_0118 + IT_0068*IT_0121 + IT_0030*IT_0124 + IT_0125);
    return create_ccomplex_return(IT_0126);
}

