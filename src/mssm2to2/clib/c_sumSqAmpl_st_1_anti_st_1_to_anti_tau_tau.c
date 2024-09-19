#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_1_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_1_to_anti_tau_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
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
    const creal_t m_st_1 = param->m_st_1;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*m_tau*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cos(alpha);
    const ccomplex_t IT_0009 = m_t*mu_h*IT_0000;
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0001*IT_0003;
    const ccomplex_t IT_0013 = IT_0008*IT_0010;
    const ccomplex_t IT_0014 = IT_0012 + IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -2);
    const ccomplex_t IT_0017 = pow(m_t, 2);
    const ccomplex_t IT_0018 = 2*U_st_00*(((-0.5)*m_W*IT_0014*(IT_0005 + -1./3
      *IT_0004*IT_0016) + IT_0000*IT_0005*IT_0008*IT_0011*IT_0017)*conj(U_st_00)
       + 1./2*IT_0003*IT_0005*IT_0009*IT_0011*conj(U_st_10)) + U_st_10*(IT_0003
      *IT_0005*IT_0009*IT_0011*conj(U_st_00) + (-1.33333333333333)*(m_W*IT_0004
      *IT_0014*IT_0016 + (-1.5)*IT_0000*IT_0005*IT_0008*IT_0011*IT_0017)*conj
      (U_st_10));
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(e_em*IT_0018 +
       1.4142135623731*A_t*IT_0008*(conj(U_st_00)*U_st_10 + U_st_00*conj(U_st_10
      )));
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = pow(m_st_1, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*e_em*m_tau*IT_0000*IT_0002
      *IT_0005*IT_0008;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = IT_0001*IT_0008;
    const ccomplex_t IT_0029 = -IT_0003*IT_0010 + IT_0028;
    const ccomplex_t IT_0030 = (-2)*U_st_00*((IT_0000*IT_0003*IT_0005*IT_0011
      *IT_0017 + 0.5*m_W*(IT_0005 + -1./3*IT_0004*IT_0016)*IT_0029)*conj(U_st_00
      ) + -1./2*IT_0005*IT_0008*IT_0009*IT_0011*conj(U_st_10)) + U_st_10*
      (IT_0005*IT_0008*IT_0009*IT_0011*conj(U_st_00) + (-1.33333333333333)*(1.5
      *IT_0000*IT_0003*IT_0005*IT_0011*IT_0017 + m_W*IT_0004*IT_0016*IT_0029)
      *conj(U_st_10));
    const ccomplex_t IT_0031 = (0 + _Complex_I*-1)*(e_em*IT_0030 + (
      -1.4142135623731)*A_t*IT_0003*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0027*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_12 + (-2)*IT_0022 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = -IT_0025 + -IT_0036;
    const ccomplex_t IT_0038 = pow(m_tau, 2);
    const ccomplex_t IT_0039 = cpow(IT_0015, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0005
      *IT_0015 + -1./3*IT_0004*IT_0039)*U_st_00*conj(U_st_00) + IT_0004*IT_0039
      *U_st_10*conj(U_st_10));
    const ccomplex_t IT_0041 = 1.33333333333333*IT_0040;
    const ccomplex_t IT_0042 = IT_0004*IT_0039;
    const ccomplex_t IT_0043 = e_em*IT_0042;
    const ccomplex_t IT_0044 = IT_0005*IT_0015;
    const ccomplex_t IT_0045 = e_em*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0043 + -IT_0045);
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0022 + -reg_prop, -1);
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
    const ccomplex_t IT_0060 = m_tau*IT_0059;
    const ccomplex_t IT_0061 = 0.666666666666667*IT_0040;
    const ccomplex_t IT_0062 = m_tau*IT_0061;
    const ccomplex_t IT_0063 = m_tau*IT_0041;
    const ccomplex_t IT_0064 = (IT_0050 + -IT_0055)*(IT_0062 + -IT_0063);
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = -(IT_0050 + -IT_0055)*(IT_0062 + -IT_0063);
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = -IT_0038;
    const ccomplex_t IT_0069 = s_34 + IT_0068;
    const ccomplex_t IT_0070 = e_em*m_tau*IT_0000*IT_0002*IT_0005*IT_0010;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = m_t*e_em*mu_h*IT_0000*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0073 = IT_0005*IT_0072;
    const ccomplex_t IT_0074 = m_t*e_em*mu_h*IT_0000*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0075 = IT_0005*IT_0074;
    const ccomplex_t IT_0076 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0077 = IT_0001*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0080 = IT_0001*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = (-0.5)*IT_0073 + 0.5*IT_0075 + (-0.5)*IT_0078 +
       0.5*IT_0081;
    const ccomplex_t IT_0083 = IT_0071*IT_0082;
    const ccomplex_t IT_0084 = cpow((-2)*s_12 + (-2)*IT_0022 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = e_em*m_tau*IT_0000*IT_0005;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = IT_0001*IT_0005*IT_0011;
    const ccomplex_t IT_0090 = IT_0072*IT_0089;
    const ccomplex_t IT_0091 = IT_0074*IT_0089;
    const ccomplex_t IT_0092 = IT_0010*IT_0076;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = IT_0010*IT_0079;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = 0.5*IT_0090 + (-0.5)*IT_0091 + (-0.5)*IT_0093 +
       0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0088*IT_0096;
    const ccomplex_t IT_0098 = pow(m_Z, 2);
    const ccomplex_t IT_0099 = cpow((-2)*s_12 + (-2)*IT_0022 + IT_0098 + 
      -reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0097*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = -IT_0086 + -IT_0101;
    const ccomplex_t IT_0103 = s_13 + s_14;
    const ccomplex_t IT_0104 = m_tau*IT_0103;
    const ccomplex_t IT_0105 = s_34 + IT_0038;
    const ccomplex_t IT_0106 = (-6)*IT_0065;
    const ccomplex_t IT_0107 = IT_0022*IT_0038;
    const ccomplex_t IT_0108 = s_14*m_tau;
    const ccomplex_t IT_0109 = (-6)*conj(IT_0052);
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = (-1.33333333333333)*IT_0110;
    const ccomplex_t IT_0113 = IT_0111*IT_0112;
    const ccomplex_t IT_0114 = cpow(s_12 + IT_0022 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = s_13*s_14;
    const ccomplex_t IT_0119 = s_34*IT_0022;
    const ccomplex_t IT_0120 = (-6)*IT_0119;
    const ccomplex_t IT_0121 = 12*IT_0118;
    const ccomplex_t IT_0122 = IT_0120 + IT_0121;
    const ccomplex_t IT_0123 = (-6)*IT_0107;
    const ccomplex_t IT_0124 = IT_0122 + IT_0123;
    const ccomplex_t IT_0125 = 6*conj(IT_0067);
    const ccomplex_t IT_0126 = IT_0117*((-12)*IT_0107 + 24*IT_0118 + (-12)
      *IT_0119) + (IT_0052 + -conj(IT_0052) + IT_0057 + -conj(IT_0057))*IT_0124 
      + -IT_0060*((-12)*IT_0037 + 12*conj(IT_0037) + 6*conj(IT_0065) + (-6)
      *IT_0067 + IT_0106 + IT_0125);
    const ccomplex_t IT_0127 = s_13*m_tau;
    const ccomplex_t IT_0128 = IT_0037*IT_0060;
    const ccomplex_t IT_0129 = (-6)*IT_0052;
    const ccomplex_t IT_0130 = 0.166666666666667*IT_0108;
    const ccomplex_t IT_0131 = 6*s_34;
    const ccomplex_t IT_0132 = IT_0037*IT_0069;
    const ccomplex_t IT_0133 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0134 = IT_0102*IT_0105;
    const ccomplex_t IT_0135 = conj(IT_0037)*(IT_0037*(12*s_34 + (-12)*IT_0038
      ) + 6*(IT_0052 + IT_0057)*IT_0060 + 6*(IT_0065 + IT_0067)*IT_0069) + conj
      (IT_0102)*((12*s_34 + 12*IT_0038)*IT_0102 + (-6)*(IT_0052 + -IT_0057)
      *IT_0104 + 6*IT_0105*(IT_0067 + 0.166666666666667*IT_0106)) + (IT_0057
      *IT_0107 + IT_0067*IT_0108)*IT_0109 + -IT_0117*IT_0126 + (-6)*conj(IT_0052
      )*(IT_0102*IT_0104 + (-0.166666666666667)*IT_0052*IT_0122 + -IT_0065
      *IT_0127 + -IT_0128) + 6*conj(IT_0057)*(IT_0102*IT_0104 +
       0.166666666666667*IT_0057*IT_0122 + IT_0067*IT_0127 + IT_0128 +
       0.166666666666667*IT_0107*IT_0129 + IT_0106*IT_0130) + 6*conj(IT_0065)*
      (IT_0052*IT_0127 + 0.166666666666667*IT_0065*IT_0131 + IT_0132) + IT_0133*
      (IT_0038*IT_0067 + IT_0057*IT_0108 + IT_0134) + IT_0125*(0.166666666666667
      *IT_0038*IT_0106 + IT_0057*IT_0127 + IT_0129*IT_0130 + 0.166666666666667
      *IT_0067*IT_0131 + IT_0132 + IT_0134);
    return create_ccomplex_return(IT_0135);
}

