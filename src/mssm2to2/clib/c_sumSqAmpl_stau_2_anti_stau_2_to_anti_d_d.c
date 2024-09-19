#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_stau_2_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_stau_2_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
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
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*U_stau_01*conj(U_stau_01) + IT_0001*IT_0002*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0005 = (-2)*IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0007 = 0.333333333333333*IT_0006;
    const ccomplex_t IT_0008 = pow(m_stau_2, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = pow(m_d, 2);
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = m_d*s_14;
    const ccomplex_t IT_0017 = -IT_0004;
    const ccomplex_t IT_0018 = m_d*IT_0017;
    const ccomplex_t IT_0019 = m_d*IT_0005;
    const ccomplex_t IT_0020 = IT_0001*IT_0002;
    const ccomplex_t IT_0021 = e_em*IT_0020;
    const ccomplex_t IT_0022 = IT_0000*IT_0003;
    const ccomplex_t IT_0023 = e_em*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0021 + 3*IT_0023);
    const ccomplex_t IT_0025 = (-0.166666666666667)*IT_0024;
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (IT_0018 + -IT_0019)*(IT_0011 + -IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0005*IT_0027;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = (-6)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_12 + (-2)*IT_0008 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0034 = pow(m_W, -1);
    const ccomplex_t IT_0035 = cos(beta);
    const ccomplex_t IT_0036 = cpow(IT_0035, -1);
    const ccomplex_t IT_0037 = cos(alpha);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_d*e_em*IT_0003*IT_0034
      *IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = mu_h*m_tau*IT_0034;
    const ccomplex_t IT_0041 = sin(alpha);
    const ccomplex_t IT_0042 = sin(beta);
    const ccomplex_t IT_0043 = -IT_0041*IT_0042;
    const ccomplex_t IT_0044 = IT_0035*IT_0037;
    const ccomplex_t IT_0045 = IT_0043 + IT_0044;
    const ccomplex_t IT_0046 = cpow(IT_0000, -2);
    const ccomplex_t IT_0047 = pow(m_tau, 2);
    const ccomplex_t IT_0048 = (-2)*U_stau_01*(((-0.5)*m_W*IT_0045*(IT_0003 + 
      -IT_0002*IT_0046) + IT_0003*IT_0034*IT_0036*IT_0037*IT_0047)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0036*IT_0040*IT_0041*conj(U_stau_11)) +
       U_stau_11*(IT_0003*IT_0036*IT_0040*IT_0041*conj(U_stau_01) + 2*(m_W
      *IT_0002*IT_0045*IT_0046 + -IT_0003*IT_0034*IT_0036*IT_0037*IT_0047)*conj
      (U_stau_11));
    const ccomplex_t IT_0049 = (0 + _Complex_I*-1)*(e_em*IT_0048 + (
      -1.4142135623731)*A_tau*IT_0037*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0039*IT_0050;
    const ccomplex_t IT_0052 = IT_0033*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*m_d*e_em*IT_0003*IT_0034
      *IT_0036*IT_0041;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = IT_0035*IT_0041;
    const ccomplex_t IT_0058 = IT_0037*IT_0042;
    const ccomplex_t IT_0059 = IT_0057 + IT_0058;
    const ccomplex_t IT_0060 = 2*U_stau_01*((IT_0003*IT_0034*IT_0036*IT_0041
      *IT_0047 + (-0.5)*m_W*(IT_0003 + -IT_0002*IT_0046)*IT_0059)*conj(U_stau_01
      ) + 1./2*IT_0003*IT_0036*IT_0037*IT_0040*conj(U_stau_11)) + U_stau_11*
      (IT_0003*IT_0036*IT_0037*IT_0040*conj(U_stau_01) + 2*(IT_0003*IT_0034
      *IT_0036*IT_0041*IT_0047 + -m_W*IT_0002*IT_0046*IT_0059)*conj(U_stau_11));
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(e_em*IT_0060 +
       1.4142135623731*A_tau*IT_0041*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = IT_0056*IT_0062;
    const ccomplex_t IT_0064 = IT_0054*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = -IT_0053 + -IT_0065;
    const ccomplex_t IT_0067 = -s_14;
    const ccomplex_t IT_0068 = s_13 + IT_0067;
    const ccomplex_t IT_0069 = m_d*IT_0068;
    const ccomplex_t IT_0070 = -(IT_0018 + -IT_0019)*(IT_0011 + -IT_0027);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = -IT_0014;
    const ccomplex_t IT_0075 = s_34 + IT_0074;
    const ccomplex_t IT_0076 = 12*s_34;
    const ccomplex_t IT_0077 = m_d*e_em*IT_0003*IT_0034*IT_0036*IT_0042;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = e_em*mu_h*m_tau*IT_0034*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0080 = IT_0003*IT_0079;
    const ccomplex_t IT_0081 = e_em*mu_h*m_tau*IT_0034*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0082 = IT_0003*IT_0081;
    const ccomplex_t IT_0083 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0084 = IT_0042*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0087 = IT_0042*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (-0.5)*IT_0080 + 0.5*IT_0082 + (-0.5)*IT_0085 +
       0.5*IT_0088;
    const ccomplex_t IT_0090 = IT_0078*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_12 + (-2)*IT_0008 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = m_d*e_em*IT_0003*IT_0034;
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = IT_0003*IT_0036*IT_0042;
    const ccomplex_t IT_0097 = IT_0079*IT_0096;
    const ccomplex_t IT_0098 = IT_0081*IT_0096;
    const ccomplex_t IT_0099 = IT_0035*IT_0083;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = IT_0035*IT_0086;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0097 + 0.5*IT_0098 + 0.5*IT_0100 + (
      -0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0095*IT_0103;
    const ccomplex_t IT_0105 = pow(m_Z, 2);
    const ccomplex_t IT_0106 = cpow((-2)*s_12 + (-2)*IT_0008 + IT_0105 + 
      -reg_prop, -1);
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = -IT_0093 + -IT_0108;
    const ccomplex_t IT_0110 = s_34 + IT_0014;
    const ccomplex_t IT_0111 = (-6)*conj(IT_0073);
    const ccomplex_t IT_0112 = s_13 + s_14;
    const ccomplex_t IT_0113 = m_d*IT_0112;
    const ccomplex_t IT_0114 = -conj(IT_0013);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0116 = 2*IT_0115;
    const ccomplex_t IT_0117 = (-0.333333333333333)*IT_0115;
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = cpow(s_12 + IT_0008 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0120 = IT_0118*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = s_13*s_14;
    const ccomplex_t IT_0124 = s_34*IT_0008;
    const ccomplex_t IT_0125 = (-6)*IT_0073;
    const ccomplex_t IT_0126 = (-6)*IT_0124;
    const ccomplex_t IT_0127 = (-6)*IT_0015;
    const ccomplex_t IT_0128 = 12*IT_0123;
    const ccomplex_t IT_0129 = IT_0126 + IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = IT_0122*((-12)*IT_0015 + 24*IT_0123 + (-12)
      *IT_0124) + 6*IT_0069*(IT_0029 + -conj(IT_0029) + 2*IT_0066 + (-2)*conj
      (IT_0066) + -conj(IT_0073) + (-0.166666666666667)*IT_0125) + (IT_0013 +
       IT_0031 + -conj(IT_0031) + IT_0114)*IT_0129;
    const ccomplex_t IT_0131 = conj(IT_0109)*IT_0110;
    const ccomplex_t IT_0132 = m_d*s_13;
    const ccomplex_t IT_0133 = IT_0126 + IT_0128;
    const ccomplex_t IT_0134 = conj(IT_0066)*IT_0069;
    const ccomplex_t IT_0135 = (-6)*conj(IT_0031);
    const ccomplex_t IT_0136 = 0.166666666666667*IT_0016;
    const ccomplex_t IT_0137 = 6*s_34;
    const ccomplex_t IT_0138 = conj(IT_0066)*IT_0075;
    const ccomplex_t IT_0139 = (conj(IT_0013)*IT_0015 + IT_0016*conj(IT_0029))
      *IT_0032 + 6*IT_0066*((conj(IT_0013) + conj(IT_0031))*IT_0069 + (conj
      (IT_0029) + conj(IT_0073))*IT_0075 + 0.166666666666667*conj(IT_0066)*((-12
      )*IT_0014 + IT_0076)) + 6*IT_0109*(0.166666666666667*(12*IT_0014 + IT_0076
      )*conj(IT_0109) + IT_0110*(conj(IT_0029) + 0.166666666666667*IT_0111) + 
      -IT_0113*(conj(IT_0031) + IT_0114)) + -IT_0122*IT_0130 + IT_0125*(conj
      (IT_0013)*IT_0016 + IT_0014*conj(IT_0029) + IT_0131) + (-6)*IT_0031*(conj
      (IT_0109)*IT_0113 + -conj(IT_0073)*IT_0132 + (-0.166666666666667)*conj
      (IT_0031)*IT_0133 + -IT_0134) + 6*IT_0013*(conj(IT_0109)*IT_0113 + conj
      (IT_0029)*IT_0132 + 0.166666666666667*conj(IT_0013)*IT_0133 + IT_0134 +
       0.166666666666667*IT_0015*IT_0135 + IT_0111*IT_0136) + 6*IT_0029*
      (0.166666666666667*IT_0014*IT_0111 + IT_0131 + conj(IT_0013)*IT_0132 +
       IT_0135*IT_0136 + 0.166666666666667*conj(IT_0029)*IT_0137 + IT_0138) + 6
      *IT_0073*(conj(IT_0031)*IT_0132 + 0.166666666666667*conj(IT_0073)*IT_0137 
      + IT_0138);
    return create_ccomplex_return(IT_0139);
}

