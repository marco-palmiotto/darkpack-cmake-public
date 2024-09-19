#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_e_anti_snu_e_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_e_anti_snu_e_to_anti_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = IT_0005 + IT_0010;
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = sin(alpha);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_W*e_em*IT_0011*(IT_0001
      *IT_0003 + -IT_0012*IT_0013);
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0005*IT_0013;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_W*e_em*IT_0011*(IT_0003
      *IT_0012 + IT_0001*IT_0013);
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0020 + -IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = pow(m_e, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0017 + IT_0034 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_23 + IT_0017 + IT_0034 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = cpow(IT_0008, -1);
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*IT_0005*IT_0045;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = IT_0004*IT_0045;
    const ccomplex_t IT_0049 = e_em*IT_0048;
    const ccomplex_t IT_0050 = IT_0005*IT_0008;
    const ccomplex_t IT_0051 = e_em*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0049 + -IT_0051);
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = IT_0047*IT_0056;
    const ccomplex_t IT_0058 = IT_0037 + IT_0044 + -IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0045;
    const ccomplex_t IT_0060 = IT_0054*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = IT_0047*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_e*U_d1*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0064*IT_0066;
    const ccomplex_t IT_0068 = IT_0035*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*m_e*U_d2*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = IT_0071*IT_0073;
    const ccomplex_t IT_0075 = IT_0042*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = -IT_0062 + IT_0069 + IT_0076;
    const ccomplex_t IT_0078 = -s_14;
    const ccomplex_t IT_0079 = s_13 + IT_0078;
    const ccomplex_t IT_0080 = m_e*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0082 = IT_0031*IT_0066;
    const ccomplex_t IT_0083 = IT_0035*IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0085 = IT_0039*IT_0073;
    const ccomplex_t IT_0086 = IT_0042*IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (-0.5)*IT_0046;
    const ccomplex_t IT_0088 = m_e*IT_0087;
    const ccomplex_t IT_0089 = m_e*IT_0047;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = IT_0088 + IT_0090;
    const ccomplex_t IT_0092 = IT_0056*IT_0091;
    const ccomplex_t IT_0093 = IT_0061*IT_0091;
    const ccomplex_t IT_0094 = m_e*IT_0037;
    const ccomplex_t IT_0095 = m_e*IT_0044;
    const ccomplex_t IT_0096 = IT_0083 + IT_0086 + -IT_0092 + IT_0093 + 
      -IT_0094 + -IT_0095;
    const ccomplex_t IT_0097 = IT_0032*IT_0064;
    const ccomplex_t IT_0098 = IT_0035*IT_0081*IT_0097;
    const ccomplex_t IT_0099 = IT_0040*IT_0071;
    const ccomplex_t IT_0100 = IT_0042*IT_0084*IT_0099;
    const ccomplex_t IT_0101 = m_e*IT_0076;
    const ccomplex_t IT_0102 = m_e*IT_0069;
    const ccomplex_t IT_0103 = IT_0092 + -IT_0093 + IT_0098 + IT_0100 + 
      -IT_0101 + -IT_0102;
    const ccomplex_t IT_0104 = -IT_0034;
    const ccomplex_t IT_0105 = s_34 + IT_0104;
    const ccomplex_t IT_0106 = m_e*s_13;
    const ccomplex_t IT_0107 = m_e*s_14;
    const ccomplex_t IT_0108 = IT_0017*IT_0034;
    const ccomplex_t IT_0109 = s_13*s_14;
    const ccomplex_t IT_0110 = s_34*IT_0017;
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0109 + IT_0111;
    const ccomplex_t IT_0113 = 2*s_34;
    const ccomplex_t IT_0114 = -conj(IT_0058);
    const ccomplex_t IT_0115 = (-2)*conj(IT_0096);
    const ccomplex_t IT_0116 = 2*IT_0029*((conj(IT_0058) + conj(IT_0077))
      *IT_0080 + 2*(conj(IT_0029) + 1./2*conj(IT_0096) + 1./2*conj(IT_0103))
      *IT_0105) + 2*IT_0058*(conj(IT_0029)*IT_0080 + conj(IT_0096)*IT_0106 + 
      -conj(IT_0103)*IT_0107 + -conj(IT_0077)*IT_0108 + 2*conj(IT_0058)*IT_0112)
       + (-2)*IT_0096*(IT_0034*conj(IT_0103) + -conj(IT_0029)*IT_0105 + conj
      (IT_0077)*IT_0107 + (-0.5)*conj(IT_0096)*IT_0113 + IT_0106*IT_0114) + 2
      *IT_0103*(conj(IT_0029)*IT_0105 + conj(IT_0077)*IT_0106 + 0.5*conj(IT_0103
      )*IT_0113 + IT_0107*IT_0114 + 0.5*IT_0034*IT_0115) + 2*IT_0077*(conj
      (IT_0029)*IT_0080 + conj(IT_0103)*IT_0106 + 2*conj(IT_0077)*IT_0112 +
       IT_0108*IT_0114 + 0.5*IT_0107*IT_0115);
    return create_ccomplex_return(IT_0116);
}

