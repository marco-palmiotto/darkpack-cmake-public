#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_mu_anti_snu_mu_to_anti_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_mu_anti_snu_mu_to_anti_mu_mu(
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
    const creal_t m_mu = param->m_mu;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = e_em*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0002, -1);
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(IT_0004 + -IT_0007);
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0005;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = m_mu*IT_0015;
    const ccomplex_t IT_0017 = -IT_0014;
    const ccomplex_t IT_0018 = m_mu*IT_0017;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = IT_0016 + IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0020*IT_0024;
    const ccomplex_t IT_0026 = pow(m_W, -1);
    const ccomplex_t IT_0027 = cos(beta);
    const ccomplex_t IT_0028 = cpow(IT_0027, -1);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_mu*IT_0005*IT_0026*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em*m_mu
      *IT_0005*IT_0026*IT_0028;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = pow(m_mu, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0010 + IT_0034 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = m_mu*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0005*IT_0026*IT_0028;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em*m_mu
      *IT_0005*IT_0026*IT_0028;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0010 + IT_0034 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = m_mu*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0049 = IT_0030*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0051 = IT_0035*IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0053 = IT_0040*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0055 = IT_0044*IT_0053*IT_0054;
    const ccomplex_t IT_0056 = IT_0021 + -IT_0025 + -IT_0038 + -IT_0047 +
       IT_0051 + IT_0055;
    const ccomplex_t IT_0057 = IT_0017*IT_0024;
    const ccomplex_t IT_0058 = IT_0037 + IT_0046 + -IT_0057;
    const ccomplex_t IT_0059 = m_mu*s_14;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0042*IT_0061;
    const ccomplex_t IT_0063 = IT_0044*IT_0054*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = IT_0032*IT_0065;
    const ccomplex_t IT_0067 = IT_0035*IT_0050*IT_0066;
    const ccomplex_t IT_0068 = IT_0052*IT_0061;
    const ccomplex_t IT_0069 = IT_0044*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = m_mu*IT_0070;
    const ccomplex_t IT_0072 = IT_0048*IT_0065;
    const ccomplex_t IT_0073 = IT_0035*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = m_mu*IT_0074;
    const ccomplex_t IT_0076 = -IT_0021 + IT_0025 + IT_0063 + IT_0067 + 
      -IT_0071 + -IT_0075;
    const ccomplex_t IT_0077 = (-2)*conj(IT_0076);
    const ccomplex_t IT_0078 = IT_0013*IT_0017;
    const ccomplex_t IT_0079 = IT_0070 + IT_0074 + -IT_0078;
    const ccomplex_t IT_0080 = sin(alpha);
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*e_em*m_mu*IT_0005*IT_0026
      *IT_0028*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = cpow(IT_0000, -2);
    const ccomplex_t IT_0084 = IT_0002*IT_0083;
    const ccomplex_t IT_0085 = IT_0005 + IT_0084;
    const ccomplex_t IT_0086 = cos(alpha);
    const ccomplex_t IT_0087 = sin(beta);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*m_W*e_em*IT_0085*(IT_0027
      *IT_0080 + IT_0086*IT_0087);
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = IT_0082*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*e_em*m_mu*IT_0005*IT_0026
      *IT_0028*IT_0086;
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_W*e_em*IT_0085*(IT_0027
      *IT_0086 + -IT_0080*IT_0087);
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = IT_0095*IT_0097;
    const ccomplex_t IT_0099 = cpow((-2)*s_12 + (-2)*IT_0010 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = -IT_0093 + -IT_0101;
    const ccomplex_t IT_0103 = -IT_0034;
    const ccomplex_t IT_0104 = s_34 + IT_0103;
    const ccomplex_t IT_0105 = 2*s_34;
    const ccomplex_t IT_0106 = m_mu*s_13;
    const ccomplex_t IT_0107 = -IT_0058;
    const ccomplex_t IT_0108 = -s_14;
    const ccomplex_t IT_0109 = s_13 + IT_0108;
    const ccomplex_t IT_0110 = m_mu*IT_0109;
    const ccomplex_t IT_0111 = 0.5*IT_0076;
    const ccomplex_t IT_0112 = IT_0010*IT_0034;
    const ccomplex_t IT_0113 = s_13*s_14;
    const ccomplex_t IT_0114 = s_34*IT_0010;
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = IT_0113 + IT_0115;
    const ccomplex_t IT_0117 = IT_0102*IT_0110;
    const ccomplex_t IT_0118 = (IT_0034*IT_0056 + IT_0058*IT_0059)*IT_0077 + (
      -2)*conj(IT_0056)*(IT_0034*IT_0076 + IT_0059*IT_0079 + -IT_0102*IT_0104 + 
      (-0.5)*IT_0056*IT_0105 + IT_0106*IT_0107) + 2*conj(IT_0102)*((IT_0058 +
       IT_0079)*IT_0110 + 2*IT_0104*(0.5*IT_0056 + IT_0102 + IT_0111)) + 2*conj
      (IT_0076)*(IT_0102*IT_0104 + IT_0079*IT_0106 + IT_0105*IT_0111) + (-2)
      *conj(IT_0058)*(IT_0059*IT_0076 + -IT_0056*IT_0106 + IT_0079*IT_0112 + (-2
      )*IT_0058*IT_0116 + -IT_0117) + (-2)*conj(IT_0079)*(IT_0056*IT_0059 + 
      -IT_0076*IT_0106 + -IT_0107*IT_0112 + (-2)*IT_0079*IT_0116 + -IT_0117);
    return create_ccomplex_return(IT_0118);
}

