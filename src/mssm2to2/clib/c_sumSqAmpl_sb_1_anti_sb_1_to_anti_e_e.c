#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_anti_e_e(
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
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_00*conj(U_sb_00) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0005 = (-0.333333333333333)*IT_0004;
    const ccomplex_t IT_0006 = m_e*IT_0005;
    const ccomplex_t IT_0007 = (-0.666666666666667)*IT_0004;
    const ccomplex_t IT_0008 = m_e*IT_0007;
    const ccomplex_t IT_0009 = pow(m_sb_1, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0001*IT_0002;
    const ccomplex_t IT_0012 = e_em*IT_0011;
    const ccomplex_t IT_0013 = IT_0000*IT_0003;
    const ccomplex_t IT_0014 = e_em*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0012 + -IT_0014);
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0020 = IT_0010*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (IT_0006 + -IT_0008)*(IT_0018 + -IT_0021);
    const ccomplex_t IT_0023 = pow(m_e, 2);
    const ccomplex_t IT_0024 = pow(m_W, -1);
    const ccomplex_t IT_0025 = cos(beta);
    const ccomplex_t IT_0026 = cpow(IT_0025, -1);
    const ccomplex_t IT_0027 = sin(beta);
    const ccomplex_t IT_0028 = m_e*e_em*IT_0003*IT_0024*IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = m_b*e_em*mu_h*IT_0024*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0031 = IT_0003*IT_0030;
    const ccomplex_t IT_0032 = m_b*e_em*mu_h*IT_0024*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = cpow(V_tb, 2);
    const ccomplex_t IT_0035 = A_b*IT_0034*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0036 = IT_0027*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = A_b*IT_0034*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0039 = IT_0027*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (-0.5)*IT_0031 + 0.5*IT_0033 + (-0.5)*IT_0037 +
       0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0029*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + (-2)*IT_0009 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = m_e*e_em*IT_0003*IT_0024;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0003*IT_0026*IT_0027;
    const ccomplex_t IT_0049 = IT_0030*IT_0048;
    const ccomplex_t IT_0050 = IT_0032*IT_0048;
    const ccomplex_t IT_0051 = IT_0025*IT_0035;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = IT_0025*IT_0038;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (-0.5)*IT_0049 + 0.5*IT_0050 + 0.5*IT_0052 + (
      -0.5)*IT_0054;
    const ccomplex_t IT_0056 = IT_0047*IT_0055;
    const ccomplex_t IT_0057 = pow(m_Z, 2);
    const ccomplex_t IT_0058 = cpow((-2)*s_12 + (-2)*IT_0009 + IT_0057 + 
      -reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0056*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = -IT_0045 + -IT_0060;
    const ccomplex_t IT_0062 = m_e*s_13;
    const ccomplex_t IT_0063 = IT_0007*IT_0021;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = sin(alpha);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*m_e*e_em*IT_0003*IT_0024
      *IT_0026*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = m_b*mu_h*IT_0024;
    const ccomplex_t IT_0069 = cos(alpha);
    const ccomplex_t IT_0070 = IT_0025*IT_0065;
    const ccomplex_t IT_0071 = IT_0027*IT_0069;
    const ccomplex_t IT_0072 = IT_0070 + IT_0071;
    const ccomplex_t IT_0073 = cpow(IT_0000, -2);
    const ccomplex_t IT_0074 = pow(m_b, 2);
    const ccomplex_t IT_0075 = 2*U_sb_00*(((-0.5)*m_W*IT_0072*(IT_0003 + 1./3
      *IT_0002*IT_0073) + IT_0003*IT_0024*IT_0026*IT_0065*IT_0074)*conj(U_sb_00)
       + 1./2*IT_0003*IT_0026*IT_0068*IT_0069*conj(U_sb_10)) + U_sb_10*(IT_0003
      *IT_0026*IT_0068*IT_0069*conj(U_sb_00) + (-0.666666666666667)*(m_W*IT_0002
      *IT_0072*IT_0073 + (-3)*IT_0003*IT_0024*IT_0026*IT_0065*IT_0074)*conj
      (U_sb_10));
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(e_em*IT_0075 +
       1.4142135623731*A_b*IT_0034*IT_0065*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0067*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*m_e*e_em*IT_0003*IT_0024
      *IT_0026*IT_0069;
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = -IT_0027*IT_0065;
    const ccomplex_t IT_0085 = IT_0025*IT_0069;
    const ccomplex_t IT_0086 = IT_0084 + IT_0085;
    const ccomplex_t IT_0087 = (-2)*U_sb_00*((IT_0003*IT_0024*IT_0026*IT_0069
      *IT_0074 + (-0.5)*m_W*(IT_0003 + 1./3*IT_0002*IT_0073)*IT_0086)*conj
      (U_sb_00) + -1./2*IT_0003*IT_0026*IT_0065*IT_0068*conj(U_sb_10)) + U_sb_10
      *(IT_0003*IT_0026*IT_0065*IT_0068*conj(U_sb_00) + (-2)*(IT_0003*IT_0024
      *IT_0026*IT_0069*IT_0074 + -1./3*m_W*IT_0002*IT_0073*IT_0086)*conj(U_sb_10
      ));
    const ccomplex_t IT_0088 = (0 + _Complex_I*-1)*(e_em*IT_0087 + (
      -1.4142135623731)*A_b*IT_0034*IT_0069*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0083*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_12 + (-2)*IT_0009 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = -IT_0081 + -IT_0093;
    const ccomplex_t IT_0095 = -IT_0023;
    const ccomplex_t IT_0096 = s_34 + IT_0095;
    const ccomplex_t IT_0097 = 6*IT_0022;
    const ccomplex_t IT_0098 = 12*s_34;
    const ccomplex_t IT_0099 = 6*conj(IT_0022);
    const ccomplex_t IT_0100 = -IT_0022;
    const ccomplex_t IT_0101 = 6*IT_0096;
    const ccomplex_t IT_0102 = IT_0007*IT_0018;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = conj(IT_0064) + conj(IT_0103);
    const ccomplex_t IT_0105 = -s_14;
    const ccomplex_t IT_0106 = s_13 + IT_0105;
    const ccomplex_t IT_0107 = m_e*IT_0106;
    const ccomplex_t IT_0108 = 6*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = 0.666666666666667*IT_0109;
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = cpow(s_12 + IT_0009 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = s_13*s_14;
    const ccomplex_t IT_0118 = IT_0009*IT_0023;
    const ccomplex_t IT_0119 = s_34*IT_0009;
    const ccomplex_t IT_0120 = (-6)*IT_0022;
    const ccomplex_t IT_0121 = IT_0116*(conj(IT_0100)*IT_0108 + -IT_0116*(24
      *IT_0117 + (-12)*IT_0118 + (-12)*IT_0119) + IT_0107*((-12)*IT_0094 + 12
      *conj(IT_0094) + IT_0099 + (-6)*IT_0100 + IT_0120));
    const ccomplex_t IT_0122 = m_e*s_14;
    const ccomplex_t IT_0123 = 6*IT_0062;
    const ccomplex_t IT_0124 = (-6)*conj(IT_0022);
    const ccomplex_t IT_0125 = (-6)*s_34 + (-6)*IT_0023;
    const ccomplex_t IT_0126 = (-6)*IT_0122;
    const ccomplex_t IT_0127 = s_13 + s_14;
    const ccomplex_t IT_0128 = m_e*IT_0127;
    const ccomplex_t IT_0129 = 6*IT_0128;
    const ccomplex_t IT_0130 = (-6)*IT_0128;
    const ccomplex_t IT_0131 = (-6)*IT_0119;
    const ccomplex_t IT_0132 = 12*IT_0117;
    const ccomplex_t IT_0133 = IT_0131 + IT_0132;
    const ccomplex_t IT_0134 = (-6)*IT_0118;
    const ccomplex_t IT_0135 = conj(IT_0094)*IT_0108;
    const ccomplex_t IT_0136 = IT_0131 + IT_0132 + IT_0134;
    const ccomplex_t IT_0137 = -IT_0116*IT_0136;
    const ccomplex_t IT_0138 = IT_0116*IT_0136;
    const ccomplex_t IT_0139 = 6*s_34*IT_0022*conj(IT_0022) + ((s_34 + IT_0023
      )*conj(IT_0061) + IT_0062*conj(IT_0064) + conj(IT_0094)*IT_0096)*IT_0097 +
       IT_0094*(conj(IT_0094)*((-12)*IT_0023 + IT_0098) + IT_0096*IT_0099 + conj
      (IT_0100)*IT_0101 + IT_0104*IT_0108) + IT_0121 + IT_0120*(IT_0023*conj
      (IT_0100) + conj(IT_0103)*IT_0122) + 6*IT_0100*(s_34*conj(IT_0100) +
       0.166666666666667*conj(IT_0094)*IT_0101 + 0.166666666666667*conj(IT_0103)
      *IT_0123 + 0.166666666666667*IT_0023*IT_0124 + 0.166666666666667*conj
      (IT_0061)*IT_0125 + 0.166666666666667*conj(IT_0064)*IT_0126) + IT_0061*
      (conj(IT_0061)*(12*IT_0023 + IT_0098) + (s_34 + IT_0023)*IT_0099 + conj
      (IT_0100)*IT_0125 + conj(IT_0064)*IT_0129 + conj(IT_0103)*IT_0130) +
       IT_0103*(conj(IT_0100)*IT_0123 + IT_0122*IT_0124 + conj(IT_0061)*IT_0130 
      + conj(IT_0103)*IT_0133 + conj(IT_0064)*IT_0134 + IT_0135 + IT_0137) +
       IT_0064*(IT_0062*IT_0099 + conj(IT_0100)*IT_0126 + conj(IT_0061)*IT_0129 
      + conj(IT_0064)*IT_0133 + conj(IT_0103)*IT_0134 + IT_0135 + IT_0137) +
       IT_0104*IT_0138;
    return create_ccomplex_return(IT_0139);
}

