#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_2_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_2_to_anti_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_e = param->m_e;
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
    const ccomplex_t IT_0005 = (-0.333333333333333)*IT_0004;
    const ccomplex_t IT_0006 = m_e*IT_0005;
    const ccomplex_t IT_0007 = (-0.666666666666667)*IT_0004;
    const ccomplex_t IT_0008 = m_e*IT_0007;
    const ccomplex_t IT_0009 = pow(m_sb_1, 2);
    const ccomplex_t IT_0010 = pow(m_sb_2, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0009 + -IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = IT_0001*IT_0002;
    const ccomplex_t IT_0016 = e_em*IT_0015;
    const ccomplex_t IT_0017 = IT_0000*IT_0003;
    const ccomplex_t IT_0018 = e_em*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0016 + -IT_0018);
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0011*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (IT_0006 + -IT_0008)*(IT_0014 + -IT_0022);
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = m_e*s_13;
    const ccomplex_t IT_0026 = IT_0007*IT_0014;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_12 + -IT_0009 + -IT_0010 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = cos(beta);
    const ccomplex_t IT_0031 = cpow(IT_0030, -1);
    const ccomplex_t IT_0032 = cos(alpha);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_e*e_em*IT_0003*IT_0029
      *IT_0031*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = cpow(V_tb, 2);
    const ccomplex_t IT_0036 = m_b*mu_h*IT_0029;
    const ccomplex_t IT_0037 = sin(alpha);
    const ccomplex_t IT_0038 = sin(beta);
    const ccomplex_t IT_0039 = -IT_0037*IT_0038;
    const ccomplex_t IT_0040 = IT_0030*IT_0032;
    const ccomplex_t IT_0041 = IT_0039 + IT_0040;
    const ccomplex_t IT_0042 = cpow(IT_0000, -2);
    const ccomplex_t IT_0043 = pow(m_b, 2);
    const ccomplex_t IT_0044 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0041*(IT_0003 + 1.
      /3*IT_0002*IT_0042) + IT_0003*IT_0029*IT_0031*IT_0032*IT_0043)*conj
      (U_sb_01) + -1./2*IT_0003*IT_0031*IT_0036*IT_0037*conj(U_sb_11)) + U_sb_10
      *(IT_0003*IT_0031*IT_0036*IT_0037*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0002*IT_0041*IT_0042 + (-3)*IT_0003*IT_0029*IT_0031*IT_0032*IT_0043)
      *conj(U_sb_11));
    const ccomplex_t IT_0045 = (0 + _Complex_I*-1)*(e_em*IT_0044 + (
      -1.4142135623731)*A_b*IT_0032*IT_0035*(conj(U_sb_01)*U_sb_10 + U_sb_00
      *conj(U_sb_11)));
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = IT_0034*IT_0046;
    const ccomplex_t IT_0048 = IT_0028*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*m_e*e_em*IT_0003*IT_0029
      *IT_0031*IT_0037;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = IT_0030*IT_0037;
    const ccomplex_t IT_0053 = IT_0032*IT_0038;
    const ccomplex_t IT_0054 = IT_0052 + IT_0053;
    const ccomplex_t IT_0055 = 2*U_sb_00*((IT_0003*IT_0029*IT_0031*IT_0037
      *IT_0043 + (-0.5)*m_W*(IT_0003 + 1./3*IT_0002*IT_0042)*IT_0054)*conj
      (U_sb_01) + 1./2*IT_0003*IT_0031*IT_0032*IT_0036*conj(U_sb_11)) + U_sb_10*
      (IT_0003*IT_0031*IT_0032*IT_0036*conj(U_sb_01) + 2*(IT_0003*IT_0029
      *IT_0031*IT_0037*IT_0043 + -1./3*m_W*IT_0002*IT_0042*IT_0054)*conj(U_sb_11
      ));
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(e_em*IT_0055 +
       1.4142135623731*A_b*IT_0035*IT_0037*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = IT_0051*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0009 + -IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = -IT_0049 + -IT_0061;
    const ccomplex_t IT_0063 = pow(m_e, 2);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = s_34 + IT_0064;
    const ccomplex_t IT_0066 = m_e*e_em*IT_0003*IT_0029*IT_0031*IT_0038;
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = m_b*e_em*mu_h*IT_0029*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0069 = IT_0003*IT_0068;
    const ccomplex_t IT_0070 = m_b*e_em*mu_h*IT_0029*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0071 = IT_0003*IT_0070;
    const ccomplex_t IT_0072 = A_b*IT_0035*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0073 = IT_0038*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = A_b*IT_0035*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0076 = IT_0038*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0069 + 0.5*IT_0071 + (-0.5)*IT_0074 +
       0.5*IT_0077;
    const ccomplex_t IT_0079 = IT_0067*IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_12 + -IT_0009 + -IT_0010 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = pow(m_Z, 2);
    const ccomplex_t IT_0084 = cpow((-2)*s_12 + -IT_0009 + -IT_0010 + IT_0083 
      + -reg_prop, -1);
    const ccomplex_t IT_0085 = m_e*e_em*IT_0003*IT_0029;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = IT_0003*IT_0031*IT_0038;
    const ccomplex_t IT_0088 = IT_0068*IT_0087;
    const ccomplex_t IT_0089 = IT_0070*IT_0087;
    const ccomplex_t IT_0090 = IT_0030*IT_0072;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = IT_0030*IT_0075;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = (-0.5)*IT_0088 + 0.5*IT_0089 + 0.5*IT_0091 + (
      -0.5)*IT_0093;
    const ccomplex_t IT_0095 = IT_0086*IT_0094;
    const ccomplex_t IT_0096 = IT_0084*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = -IT_0082 + -IT_0097;
    const ccomplex_t IT_0099 = s_34 + IT_0063;
    const ccomplex_t IT_0100 = 6*IT_0024;
    const ccomplex_t IT_0101 = 12*s_34;
    const ccomplex_t IT_0102 = -IT_0024;
    const ccomplex_t IT_0103 = 6*conj(IT_0024);
    const ccomplex_t IT_0104 = IT_0007*IT_0022;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = -s_14;
    const ccomplex_t IT_0107 = s_13 + IT_0106;
    const ccomplex_t IT_0108 = m_e*IT_0107;
    const ccomplex_t IT_0109 = 6*IT_0108;
    const ccomplex_t IT_0110 = s_13 + s_14;
    const ccomplex_t IT_0111 = m_e*IT_0110;
    const ccomplex_t IT_0112 = 6*IT_0111;
    const ccomplex_t IT_0113 = (-6)*IT_0111;
    const ccomplex_t IT_0114 = s_34*IT_0009;
    const ccomplex_t IT_0115 = (-6)*IT_0114;
    const ccomplex_t IT_0116 = s_13*s_14;
    const ccomplex_t IT_0117 = 12*IT_0116;
    const ccomplex_t IT_0118 = IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = conj(IT_0062)*IT_0109;
    const ccomplex_t IT_0120 = m_e*s_14;
    const ccomplex_t IT_0121 = (-6)*IT_0120;
    const ccomplex_t IT_0122 = IT_0009*IT_0063;
    const ccomplex_t IT_0123 = (-6)*IT_0122;
    const ccomplex_t IT_0124 = (-6)*IT_0024;
    const ccomplex_t IT_0125 = 6*IT_0025;
    const ccomplex_t IT_0126 = (-6)*conj(IT_0024);
    const ccomplex_t IT_0127 = 6*s_34*IT_0024*conj(IT_0024) + (IT_0025*conj
      (IT_0027) + conj(IT_0062)*IT_0065 + conj(IT_0098)*IT_0099)*IT_0100 +
       IT_0062*(conj(IT_0062)*((-12)*IT_0063 + IT_0101) + 6*IT_0065*(conj
      (IT_0102) + 0.166666666666667*IT_0103) + (conj(IT_0027) + conj(IT_0105))
      *IT_0109) + IT_0098*(conj(IT_0098)*(12*IT_0063 + IT_0101) + (-6)*IT_0099*
      (conj(IT_0102) + (-0.166666666666667)*IT_0103) + conj(IT_0027)*IT_0112 +
       conj(IT_0105)*IT_0113) + IT_0027*(IT_0025*IT_0103 + conj(IT_0098)*IT_0112
       + conj(IT_0027)*IT_0118 + IT_0119 + conj(IT_0102)*IT_0121 + conj(IT_0105)
      *IT_0123) + (IT_0063*conj(IT_0102) + conj(IT_0105)*IT_0120)*IT_0124 + 6
      *IT_0102*(conj(IT_0062)*IT_0065 + -conj(IT_0098)*IT_0099 + s_34*conj
      (IT_0102) + 0.166666666666667*conj(IT_0027)*IT_0121 + 0.166666666666667
      *conj(IT_0105)*IT_0125 + 0.166666666666667*IT_0063*IT_0126) + IT_0105*
      (conj(IT_0098)*IT_0113 + conj(IT_0105)*IT_0118 + IT_0119 + conj(IT_0027)
      *IT_0123 + conj(IT_0102)*IT_0125 + IT_0120*IT_0126);
    return create_ccomplex_return(IT_0127);
}

