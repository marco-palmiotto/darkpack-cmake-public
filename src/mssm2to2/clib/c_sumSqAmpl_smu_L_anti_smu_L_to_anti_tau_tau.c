#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_smu_L_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_smu_L_to_anti_tau_tau(
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
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003);
    const ccomplex_t IT_0005 = (-0.5)*IT_0004;
    const ccomplex_t IT_0006 = m_tau*IT_0005;
    const ccomplex_t IT_0007 = -IT_0004;
    const ccomplex_t IT_0008 = m_tau*IT_0007;
    const ccomplex_t IT_0009 = IT_0001*IT_0002;
    const ccomplex_t IT_0010 = e_em*IT_0009;
    const ccomplex_t IT_0011 = IT_0000*IT_0003;
    const ccomplex_t IT_0012 = e_em*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0010 + -IT_0012);
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = pow(m_smu_L, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (IT_0006 + -IT_0008)*(IT_0018 + -IT_0021);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = IT_0007*IT_0018;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0007*IT_0021;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = s_13*m_tau;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = cos(beta);
    const ccomplex_t IT_0031 = cpow(IT_0030, -1);
    const ccomplex_t IT_0032 = sin(alpha);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*e_em*m_tau*IT_0003*IT_0029
      *IT_0031*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = pow(m_mu, 2);
    const ccomplex_t IT_0036 = cos(alpha);
    const ccomplex_t IT_0037 = sin(beta);
    const ccomplex_t IT_0038 = cpow(IT_0000, -2);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0029
      *IT_0031*IT_0032*IT_0035 + (-0.5)*m_W*(IT_0030*IT_0032 + IT_0036*IT_0037)*
      (IT_0003 + -IT_0002*IT_0038));
    const ccomplex_t IT_0040 = IT_0034*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*e_em*m_tau*IT_0003*IT_0029
      *IT_0031*IT_0036;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0029
      *IT_0031*IT_0035*IT_0036 + (-0.5)*m_W*(IT_0030*IT_0036 + -IT_0032*IT_0037)
      *(IT_0003 + -IT_0002*IT_0038));
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + (-2)*IT_0015 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = -IT_0043 + -IT_0051;
    const ccomplex_t IT_0053 = pow(m_tau, 2);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = s_34 + IT_0054;
    const ccomplex_t IT_0056 = s_13*s_14;
    const ccomplex_t IT_0057 = s_34*IT_0015;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = 2*s_34;
    const ccomplex_t IT_0061 = IT_0015*IT_0053;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0059 + IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = 2*IT_0064;
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = cpow(s_12 + IT_0015 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = cpow(IT_0071, 2);
    const ccomplex_t IT_0073 = conj(IT_0025) + conj(IT_0027);
    const ccomplex_t IT_0074 = 4*IT_0071;
    const ccomplex_t IT_0075 = IT_0025 + IT_0027;
    const ccomplex_t IT_0076 = (-4)*IT_0071;
    const ccomplex_t IT_0077 = -s_14;
    const ccomplex_t IT_0078 = s_13 + IT_0077;
    const ccomplex_t IT_0079 = m_tau*IT_0078;
    const ccomplex_t IT_0080 = (-2)*IT_0023;
    const ccomplex_t IT_0081 = (-2)*conj(IT_0023);
    const ccomplex_t IT_0082 = (-2)*IT_0025;
    const ccomplex_t IT_0083 = (-2)*conj(IT_0025);
    const ccomplex_t IT_0084 = s_14*m_tau;
    const ccomplex_t IT_0085 = 2*(conj(IT_0023)*IT_0025 + IT_0023*conj(IT_0025
      ) + conj(IT_0022)*IT_0027 + IT_0022*conj(IT_0027))*IT_0028 + 4*(IT_0052*
      (0.5*conj(IT_0022) + 0.5*conj(IT_0023) + conj(IT_0052)) + 1./2*(IT_0022 +
       IT_0023)*conj(IT_0052))*IT_0055 + 4*(IT_0025*conj(IT_0025) + IT_0027*conj
      (IT_0027))*IT_0059 + (IT_0022*conj(IT_0022) + IT_0023*conj(IT_0023))
      *IT_0060 + (-8)*IT_0063*(IT_0072 + (-0.125)*IT_0073*IT_0074 + (-0.125)
      *IT_0075*IT_0076) + 2*IT_0079*(IT_0052*IT_0073 + conj(IT_0052)*IT_0075 + (
      -2)*IT_0071*(0.5*IT_0022 + (-0.5)*conj(IT_0022) + (-0.5)*conj(IT_0023) +
       IT_0052 + -conj(IT_0052) + (-0.25)*IT_0080)) + IT_0053*(conj(IT_0022)
      *IT_0080 + IT_0022*IT_0081) + IT_0061*(conj(IT_0027)*IT_0082 + IT_0027
      *IT_0083) + (conj(IT_0027)*IT_0080 + IT_0027*IT_0081 + conj(IT_0022)
      *IT_0082 + IT_0022*IT_0083)*IT_0084;
    return create_ccomplex_return(IT_0085);
}

