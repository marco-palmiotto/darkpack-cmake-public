#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_2_to_anti_u_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_2_to_anti_u_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
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
    const ccomplex_t IT_0007 = (-0.666666666666667)*IT_0006;
    const ccomplex_t IT_0008 = pow(m_stau_1, 2);
    const ccomplex_t IT_0009 = pow(m_stau_2, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0008 + -IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = IT_0005*IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = pow(m_u, 2);
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = m_u*s_14;
    const ccomplex_t IT_0018 = -IT_0004;
    const ccomplex_t IT_0019 = m_u*IT_0018;
    const ccomplex_t IT_0020 = m_u*IT_0005;
    const ccomplex_t IT_0021 = IT_0001*IT_0002;
    const ccomplex_t IT_0022 = e_em*IT_0021;
    const ccomplex_t IT_0023 = IT_0000*IT_0003;
    const ccomplex_t IT_0024 = e_em*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0022 + (-3)*IT_0024);
    const ccomplex_t IT_0026 = (-0.166666666666667)*IT_0025;
    const ccomplex_t IT_0027 = IT_0010*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = (IT_0019 + -IT_0020)*(IT_0012 + -IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = IT_0005*IT_0028;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = (-6)*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0035 = pow(m_W, -1);
    const ccomplex_t IT_0036 = e_em*mu_h*m_tau*IT_0035*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = e_em*mu_h*m_tau*IT_0035*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0039 = IT_0003*IT_0038;
    const ccomplex_t IT_0040 = sin(beta);
    const ccomplex_t IT_0041 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0045 = IT_0040*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (-0.5)*IT_0037 + 0.5*IT_0039 + (-0.5)*IT_0043 +
       0.5*IT_0046;
    const ccomplex_t IT_0048 = cos(beta);
    const ccomplex_t IT_0049 = cpow(IT_0040, -1);
    const ccomplex_t IT_0050 = m_u*e_em*IT_0003*IT_0035*IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0047*IT_0051;
    const ccomplex_t IT_0053 = IT_0034*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = m_u*e_em*IT_0003*IT_0035;
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = cpow(IT_0048, -1);
    const ccomplex_t IT_0058 = IT_0003*IT_0040*IT_0057;
    const ccomplex_t IT_0059 = IT_0036*IT_0058;
    const ccomplex_t IT_0060 = IT_0038*IT_0058;
    const ccomplex_t IT_0061 = IT_0041*IT_0048;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = IT_0044*IT_0048;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (-0.5)*IT_0059 + 0.5*IT_0060 + 0.5*IT_0062 + (
      -0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0056*IT_0065;
    const ccomplex_t IT_0067 = pow(m_Z, 2);
    const ccomplex_t IT_0068 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + IT_0067 
      + -reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = -IT_0054 + -IT_0070;
    const ccomplex_t IT_0072 = s_13 + s_14;
    const ccomplex_t IT_0073 = m_u*IT_0072;
    const ccomplex_t IT_0074 = s_34 + IT_0015;
    const ccomplex_t IT_0075 = -(IT_0019 + -IT_0020)*(IT_0012 + -IT_0028);
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0078 = 12*s_34;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0008 + -IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0080 = cos(alpha);
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_u*e_em*IT_0003*IT_0035
      *IT_0049*IT_0080;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = sin(alpha);
    const ccomplex_t IT_0084 = mu_h*m_tau*IT_0035;
    const ccomplex_t IT_0085 = IT_0048*IT_0083;
    const ccomplex_t IT_0086 = IT_0040*IT_0080;
    const ccomplex_t IT_0087 = IT_0085 + IT_0086;
    const ccomplex_t IT_0088 = cpow(IT_0000, -2);
    const ccomplex_t IT_0089 = pow(m_tau, 2);
    const ccomplex_t IT_0090 = 2*U_stau_00*(((-0.5)*m_W*IT_0087*(IT_0003 + 
      -IT_0002*IT_0088) + IT_0003*IT_0035*IT_0057*IT_0083*IT_0089)*conj
      (U_stau_01) + 1./2*IT_0003*IT_0057*IT_0080*IT_0084*conj(U_stau_11)) +
       U_stau_10*(IT_0003*IT_0057*IT_0080*IT_0084*conj(U_stau_01) + (-2)*(m_W
      *IT_0002*IT_0087*IT_0088 + -IT_0003*IT_0035*IT_0057*IT_0083*IT_0089)*conj
      (U_stau_11));
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(e_em*IT_0090 +
       1.4142135623731*A_tau*IT_0083*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0082*IT_0092;
    const ccomplex_t IT_0094 = IT_0079*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_12 + -IT_0008 + -IT_0009 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*m_u*e_em*IT_0003*IT_0035
      *IT_0049*IT_0083;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = -IT_0040*IT_0083;
    const ccomplex_t IT_0100 = IT_0048*IT_0080;
    const ccomplex_t IT_0101 = IT_0099 + IT_0100;
    const ccomplex_t IT_0102 = (-2)*U_stau_00*((IT_0003*IT_0035*IT_0057
      *IT_0080*IT_0089 + (-0.5)*m_W*(IT_0003 + -IT_0002*IT_0088)*IT_0101)*conj
      (U_stau_01) + -1./2*IT_0003*IT_0057*IT_0083*IT_0084*conj(U_stau_11)) +
       U_stau_10*(IT_0003*IT_0057*IT_0083*IT_0084*conj(U_stau_01) + (-2)*
      (IT_0003*IT_0035*IT_0057*IT_0080*IT_0089 + -m_W*IT_0002*IT_0088*IT_0101)
      *conj(U_stau_11));
    const ccomplex_t IT_0103 = (0 + _Complex_I*-1)*(e_em*IT_0102 + (
      -1.4142135623731)*A_tau*IT_0080*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = IT_0098*IT_0104;
    const ccomplex_t IT_0106 = IT_0096*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = -IT_0095 + -IT_0107;
    const ccomplex_t IT_0109 = -s_14;
    const ccomplex_t IT_0110 = s_13 + IT_0109;
    const ccomplex_t IT_0111 = m_u*IT_0110;
    const ccomplex_t IT_0112 = -IT_0015;
    const ccomplex_t IT_0113 = s_34 + IT_0112;
    const ccomplex_t IT_0114 = (-6)*IT_0076;
    const ccomplex_t IT_0115 = conj(IT_0071)*IT_0074;
    const ccomplex_t IT_0116 = m_u*s_13;
    const ccomplex_t IT_0117 = s_34*IT_0008;
    const ccomplex_t IT_0118 = (-6)*IT_0117;
    const ccomplex_t IT_0119 = s_13*s_14;
    const ccomplex_t IT_0120 = 12*IT_0119;
    const ccomplex_t IT_0121 = IT_0118 + IT_0120;
    const ccomplex_t IT_0122 = conj(IT_0108)*IT_0111;
    const ccomplex_t IT_0123 = (-6)*conj(IT_0032);
    const ccomplex_t IT_0124 = 0.166666666666667*IT_0017;
    const ccomplex_t IT_0125 = 6*s_34;
    const ccomplex_t IT_0126 = conj(IT_0108)*IT_0113;
    const ccomplex_t IT_0127 = (conj(IT_0014)*IT_0016 + IT_0017*conj(IT_0030))
      *IT_0033 + 6*IT_0071*((conj(IT_0014) + -conj(IT_0032))*IT_0073 + IT_0074*
      (conj(IT_0030) + 0.166666666666667*IT_0077) + 0.166666666666667*conj
      (IT_0071)*(12*IT_0015 + IT_0078)) + 6*IT_0108*(0.166666666666667*((-12)
      *IT_0015 + IT_0078)*conj(IT_0108) + (conj(IT_0014) + conj(IT_0032))
      *IT_0111 + (conj(IT_0030) + conj(IT_0076))*IT_0113) + IT_0114*(conj
      (IT_0014)*IT_0017 + IT_0015*conj(IT_0030) + IT_0115) + (-6)*IT_0032*(conj
      (IT_0071)*IT_0073 + -conj(IT_0076)*IT_0116 + (-0.166666666666667)*conj
      (IT_0032)*IT_0121 + -IT_0122) + 6*IT_0014*(conj(IT_0071)*IT_0073 + conj
      (IT_0030)*IT_0116 + 0.166666666666667*conj(IT_0014)*IT_0121 + IT_0122 +
       0.166666666666667*IT_0016*IT_0123 + IT_0077*IT_0124) + 6*IT_0030*
      (0.166666666666667*IT_0015*IT_0077 + IT_0115 + conj(IT_0014)*IT_0116 +
       IT_0123*IT_0124 + 0.166666666666667*conj(IT_0030)*IT_0125 + IT_0126) + 6
      *IT_0076*(conj(IT_0032)*IT_0116 + 0.166666666666667*conj(IT_0076)*IT_0125 
      + IT_0126);
    return create_ccomplex_return(IT_0127);
}

