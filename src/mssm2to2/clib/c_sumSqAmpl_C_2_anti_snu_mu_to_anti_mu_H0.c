#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_snu_mu_to_anti_mu_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_snu_mu_to_anti_mu_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_H0, 2);
    const ccomplex_t IT_0001 = m_mu*m_C_2*IT_0000;
    const ccomplex_t IT_0002 = pow(m_mu, 2);
    const ccomplex_t IT_0003 = pow(m_C_2, 2);
    const ccomplex_t IT_0004 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0005 = cpow((-2)*s_12 + IT_0002 + -IT_0003 + -IT_0004 
      + -reg_prop, -1);
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = cos(alpha);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*e_em*m_mu*IT_0006*IT_0008
      *IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_mu*IT_0006*IT_0008*IT_0011;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = IT_0005*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0002 + IT_0004 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0011;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0024 = IT_0009*IT_0023;
    const ccomplex_t IT_0025 = sin(alpha);
    const ccomplex_t IT_0026 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = IT_0011*(IT_0024 + IT_0027);
    const ccomplex_t IT_0029 = 1.4142135623731*e_em*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = IT_0022*IT_0031;
    const ccomplex_t IT_0033 = IT_0020*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0002 + IT_0004 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0036 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0037 = IT_0009*IT_0036;
    const ccomplex_t IT_0038 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0039 = IT_0025*IT_0038;
    const ccomplex_t IT_0040 = IT_0011*(IT_0037 + IT_0039);
    const ccomplex_t IT_0041 = 1.4142135623731*e_em*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0011;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = IT_0035*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = -IT_0034 + -IT_0048;
    const ccomplex_t IT_0050 = conj(IT_0019) + conj(IT_0049);
    const ccomplex_t IT_0051 = U_d1*V_Wp2;
    const ccomplex_t IT_0052 = IT_0009*IT_0051;
    const ccomplex_t IT_0053 = V_u2*U_Wm1;
    const ccomplex_t IT_0054 = IT_0025*IT_0053;
    const ccomplex_t IT_0055 = IT_0011*(IT_0052 + IT_0054);
    const ccomplex_t IT_0056 = 1.4142135623731*e_em*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0006*IT_0008*IT_0011;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = IT_0058*IT_0060;
    const ccomplex_t IT_0062 = IT_0020*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = U_d2*V_Wp2;
    const ccomplex_t IT_0065 = IT_0009*IT_0064;
    const ccomplex_t IT_0066 = V_u2*U_Wm2;
    const ccomplex_t IT_0067 = IT_0025*IT_0066;
    const ccomplex_t IT_0068 = IT_0011*(IT_0065 + IT_0067);
    const ccomplex_t IT_0069 = 1.4142135623731*e_em*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = IT_0015*IT_0071;
    const ccomplex_t IT_0073 = IT_0035*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0063 + IT_0074;
    const ccomplex_t IT_0076 = IT_0013*IT_0045;
    const ccomplex_t IT_0077 = IT_0005*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0075 + IT_0078;
    const ccomplex_t IT_0080 = IT_0019 + IT_0049;
    const ccomplex_t IT_0081 = conj(IT_0075) + conj(IT_0078);
    const ccomplex_t IT_0082 = s_34*m_C_2;
    const ccomplex_t IT_0083 = m_C_2*IT_0075;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_13 + IT_0002 + IT_0003 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0086 = cos(theta_W);
    const ccomplex_t IT_0087 = cpow(IT_0086, -2);
    const ccomplex_t IT_0088 = IT_0010*IT_0087;
    const ccomplex_t IT_0089 = IT_0011 + IT_0088;
    const ccomplex_t IT_0090 = sin(beta);
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_W*e_em*IT_0089*(IT_0007
      *IT_0009 + -IT_0025*IT_0090);
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = IT_0015*IT_0092;
    const ccomplex_t IT_0094 = IT_0085*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0097 = IT_0031*IT_0060;
    const ccomplex_t IT_0098 = IT_0020*IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0100 = IT_0015*IT_0043;
    const ccomplex_t IT_0101 = IT_0035*IT_0099*IT_0100;
    const ccomplex_t IT_0102 = m_mu*IT_0018;
    const ccomplex_t IT_0103 = -IT_0095 + IT_0098 + IT_0101 + 2*IT_0102;
    const ccomplex_t IT_0104 = (-2)*conj(IT_0103);
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = conj(IT_0084) + IT_0105;
    const ccomplex_t IT_0107 = (-2)*IT_0103;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = IT_0084 + IT_0108;
    const ccomplex_t IT_0110 = IT_0034 + IT_0048;
    const ccomplex_t IT_0111 = m_C_2*IT_0110;
    const ccomplex_t IT_0112 = IT_0045*IT_0092;
    const ccomplex_t IT_0113 = IT_0085*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0022*IT_0058;
    const ccomplex_t IT_0116 = IT_0020*IT_0096*IT_0115;
    const ccomplex_t IT_0117 = IT_0045*IT_0071;
    const ccomplex_t IT_0118 = IT_0035*IT_0099*IT_0117;
    const ccomplex_t IT_0119 = m_mu*IT_0078;
    const ccomplex_t IT_0120 = IT_0114 + -IT_0116 + -IT_0118 + (-2)*IT_0119;
    const ccomplex_t IT_0121 = (-2)*IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = IT_0111 + IT_0122;
    const ccomplex_t IT_0124 = (-2)*conj(IT_0120);
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = conj(IT_0111) + IT_0125;
    const ccomplex_t IT_0127 = m_mu*s_14;
    const ccomplex_t IT_0128 = 2*s_13;
    const ccomplex_t IT_0129 = (-2)*s_13;
    const ccomplex_t IT_0130 = m_mu*m_C_2;
    const ccomplex_t IT_0131 = s_13*IT_0000;
    const ccomplex_t IT_0132 = s_14*s_34;
    const ccomplex_t IT_0133 = (-2)*IT_0132;
    const ccomplex_t IT_0134 = IT_0131 + IT_0133;
    const ccomplex_t IT_0135 = (-2)*IT_0075;
    const ccomplex_t IT_0136 = (-2)*IT_0049;
    const ccomplex_t IT_0137 = 2*IT_0001*(IT_0050*IT_0079 + IT_0080*IT_0081) +
       2*IT_0082*(IT_0079*IT_0106 + IT_0081*IT_0109 + IT_0050*IT_0123 + IT_0080
      *IT_0126) + 2*(IT_0080*IT_0106 + IT_0050*IT_0109 + IT_0081*IT_0123 +
       IT_0079*IT_0126)*IT_0127 + (IT_0084*conj(IT_0084) + IT_0103*conj(IT_0103)
       + IT_0111*conj(IT_0111) + IT_0120*conj(IT_0120))*IT_0128 + (conj(IT_0084)
      *IT_0103 + IT_0084*conj(IT_0103) + conj(IT_0111)*IT_0120 + IT_0111*conj
      (IT_0120))*IT_0129 + 2*(0.5*IT_0104*IT_0111 + 0.5*IT_0107*conj(IT_0111) +
       conj(IT_0103)*IT_0120 + IT_0103*conj(IT_0120) + conj(IT_0084)*IT_0123 +
       IT_0084*IT_0126)*IT_0130 + (-2)*IT_0134*(IT_0081*(IT_0078 + (-0.5)
      *IT_0135) + IT_0050*(IT_0019 + (-0.5)*IT_0136));
    return create_ccomplex_return(IT_0137);
}

