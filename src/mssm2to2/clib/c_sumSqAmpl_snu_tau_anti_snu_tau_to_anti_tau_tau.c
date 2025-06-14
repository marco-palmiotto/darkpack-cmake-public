#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_tau_anti_snu_tau_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_tau_anti_snu_tau_to_anti_tau_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
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
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
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
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*m_tau*IT_0000*IT_0002
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
    const ccomplex_t IT_0017 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*m_tau*IT_0000*IT_0002
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
    const ccomplex_t IT_0030 = cpow(IT_0008, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*e_em*IT_0005*IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = IT_0004*IT_0030;
    const ccomplex_t IT_0034 = e_em*IT_0033;
    const ccomplex_t IT_0035 = IT_0005*IT_0008;
    const ccomplex_t IT_0036 = e_em*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0034 + -IT_0036);
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0032*IT_0041;
    const ccomplex_t IT_0043 = pow(m_tau, 2);
    const ccomplex_t IT_0044 = cpow(s_23 + (-0.5)*IT_0017 + (-0.5)*IT_0043 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = IT_0044*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow(s_23 + (-0.5)*IT_0017 + (-0.5)*IT_0043 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = IT_0051*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = -IT_0042 + (-0.5)*IT_0050 + (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = IT_0044*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em
      *m_tau*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = IT_0051*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0030;
    const ccomplex_t IT_0074 = IT_0039*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0032*IT_0075;
    const ccomplex_t IT_0077 = (-0.5)*IT_0065 + (-0.5)*IT_0072 + -IT_0076;
    const ccomplex_t IT_0078 = -s_14;
    const ccomplex_t IT_0079 = s_13 + IT_0078;
    const ccomplex_t IT_0080 = m_tau*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0082 = IT_0053*IT_0067;
    const ccomplex_t IT_0083 = IT_0051*IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0085 = IT_0047*IT_0060;
    const ccomplex_t IT_0086 = IT_0044*IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (-0.5)*IT_0031;
    const ccomplex_t IT_0088 = m_tau*IT_0087;
    const ccomplex_t IT_0089 = m_tau*IT_0032;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = IT_0088 + IT_0090;
    const ccomplex_t IT_0092 = IT_0041*IT_0091;
    const ccomplex_t IT_0093 = m_tau*IT_0057;
    const ccomplex_t IT_0094 = IT_0075*IT_0091;
    const ccomplex_t IT_0095 = m_tau*IT_0050;
    const ccomplex_t IT_0096 = (-0.5)*IT_0083 + (-0.5)*IT_0086 + -IT_0092 +
       0.5*IT_0093 + IT_0094 + 0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0054*IT_0069;
    const ccomplex_t IT_0098 = IT_0051*IT_0081*IT_0097;
    const ccomplex_t IT_0099 = IT_0045*IT_0062;
    const ccomplex_t IT_0100 = IT_0044*IT_0084*IT_0099;
    const ccomplex_t IT_0101 = m_tau*IT_0065;
    const ccomplex_t IT_0102 = m_tau*IT_0072;
    const ccomplex_t IT_0103 = IT_0092 + -IT_0094 + (-0.5)*IT_0098 + (-0.5)
      *IT_0100 + 0.5*IT_0101 + 0.5*IT_0102;
    const ccomplex_t IT_0104 = -IT_0043;
    const ccomplex_t IT_0105 = s_34 + IT_0104;
    const ccomplex_t IT_0106 = s_13*m_tau;
    const ccomplex_t IT_0107 = 2*s_34;
    const ccomplex_t IT_0108 = s_14*m_tau;
    const ccomplex_t IT_0109 = (-2)*conj(IT_0096);
    const ccomplex_t IT_0110 = -IT_0077;
    const ccomplex_t IT_0111 = IT_0017*IT_0043;
    const ccomplex_t IT_0112 = s_13*s_14;
    const ccomplex_t IT_0113 = s_34*IT_0017;
    const ccomplex_t IT_0114 = (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = IT_0112 + IT_0114;
    const ccomplex_t IT_0116 = IT_0029*IT_0080;
    const ccomplex_t IT_0117 = 2*conj(IT_0029)*((IT_0058 + IT_0077)*IT_0080 +
       2*(IT_0029 + 1./2*IT_0096 + 1./2*IT_0103)*IT_0105) + 2*conj(IT_0096)*
      (IT_0029*IT_0105 + IT_0058*IT_0106 + 0.5*IT_0096*IT_0107) + (IT_0043
      *IT_0103 + IT_0077*IT_0108)*IT_0109 + (-2)*conj(IT_0103)*(IT_0043*IT_0096 
      + -IT_0029*IT_0105 + (-0.5)*IT_0103*IT_0107 + IT_0058*IT_0108 + IT_0106
      *IT_0110) + 2*conj(IT_0058)*(IT_0096*IT_0106 + -IT_0103*IT_0108 + IT_0110
      *IT_0111 + 2*IT_0058*IT_0115 + IT_0116) + 2*conj(IT_0077)*(IT_0103*IT_0106
       + -IT_0096*IT_0108 + -IT_0058*IT_0111 + 2*IT_0077*IT_0115 + IT_0116);
    return create_ccomplex_return(IT_0117);
}

