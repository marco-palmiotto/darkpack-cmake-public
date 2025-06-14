#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_se_R_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_se_R_to_anti_tau_tau(
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
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0004 = -IT_0003;
    const ccomplex_t IT_0005 = m_tau*IT_0004;
    const ccomplex_t IT_0006 = (-2)*IT_0003;
    const ccomplex_t IT_0007 = m_tau*IT_0006;
    const ccomplex_t IT_0008 = IT_0001*IT_0002;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0002, -1);
    const ccomplex_t IT_0011 = IT_0000*IT_0010;
    const ccomplex_t IT_0012 = e_em*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0009 + -IT_0012);
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = pow(m_se_R, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0003*IT_0016;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (IT_0005 + -IT_0007)*(IT_0018 + -IT_0020);
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0006*IT_0018;
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = IT_0006*IT_0020;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = s_13*m_tau;
    const ccomplex_t IT_0028 = pow(m_W, -1);
    const ccomplex_t IT_0029 = cos(beta);
    const ccomplex_t IT_0030 = cpow(IT_0029, -1);
    const ccomplex_t IT_0031 = sin(alpha);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*m_tau*IT_0010*IT_0028
      *IT_0030*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = cpow(IT_0000, -2);
    const ccomplex_t IT_0035 = cos(alpha);
    const ccomplex_t IT_0036 = sin(beta);
    const ccomplex_t IT_0037 = pow(m_e, 2);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*e_em*(m_W*IT_0002*IT_0034*
      (IT_0029*IT_0031 + IT_0035*IT_0036) + -IT_0010*IT_0028*IT_0030*IT_0031
      *IT_0037);
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = IT_0033*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*e_em*m_tau*IT_0010*IT_0028
      *IT_0030*IT_0035;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*(m_W*IT_0002*IT_0034*
      (IT_0029*IT_0035 + -IT_0031*IT_0036) + -IT_0010*IT_0028*IT_0030*IT_0035
      *IT_0037);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + (-2)*IT_0015 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = -IT_0043 + -IT_0050;
    const ccomplex_t IT_0052 = pow(m_tau, 2);
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = s_34 + IT_0053;
    const ccomplex_t IT_0055 = s_13*s_14;
    const ccomplex_t IT_0056 = s_34*IT_0015;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0055 + IT_0057;
    const ccomplex_t IT_0059 = 2*s_34;
    const ccomplex_t IT_0060 = IT_0015*IT_0052;
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0058 + IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = 2*IT_0063;
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = cpow(s_12 + IT_0015 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = cpow(IT_0070, 2);
    const ccomplex_t IT_0072 = conj(IT_0024) + conj(IT_0026);
    const ccomplex_t IT_0073 = 4*IT_0070;
    const ccomplex_t IT_0074 = IT_0024 + IT_0026;
    const ccomplex_t IT_0075 = (-4)*IT_0070;
    const ccomplex_t IT_0076 = -s_14;
    const ccomplex_t IT_0077 = s_13 + IT_0076;
    const ccomplex_t IT_0078 = m_tau*IT_0077;
    const ccomplex_t IT_0079 = (-2)*IT_0022;
    const ccomplex_t IT_0080 = (-2)*conj(IT_0022);
    const ccomplex_t IT_0081 = (-2)*IT_0024;
    const ccomplex_t IT_0082 = (-2)*conj(IT_0024);
    const ccomplex_t IT_0083 = s_14*m_tau;
    const ccomplex_t IT_0084 = 2*(conj(IT_0022)*IT_0024 + IT_0022*conj(IT_0024
      ) + conj(IT_0021)*IT_0026 + IT_0021*conj(IT_0026))*IT_0027 + 4*(IT_0051*
      (0.5*conj(IT_0021) + 0.5*conj(IT_0022) + conj(IT_0051)) + 1./2*(IT_0021 +
       IT_0022)*conj(IT_0051))*IT_0054 + 4*(IT_0024*conj(IT_0024) + IT_0026*conj
      (IT_0026))*IT_0058 + (IT_0021*conj(IT_0021) + IT_0022*conj(IT_0022))
      *IT_0059 + (-8)*IT_0062*(IT_0071 + (-0.125)*IT_0072*IT_0073 + (-0.125)
      *IT_0074*IT_0075) + 2*IT_0078*(IT_0051*IT_0072 + conj(IT_0051)*IT_0074 + (
      -2)*IT_0070*(0.5*IT_0021 + (-0.5)*conj(IT_0021) + (-0.5)*conj(IT_0022) +
       IT_0051 + -conj(IT_0051) + (-0.25)*IT_0079)) + IT_0052*(conj(IT_0021)
      *IT_0079 + IT_0021*IT_0080) + IT_0060*(conj(IT_0026)*IT_0081 + IT_0026
      *IT_0082) + (conj(IT_0026)*IT_0079 + IT_0026*IT_0080 + conj(IT_0021)
      *IT_0081 + IT_0021*IT_0082)*IT_0083;
    return create_ccomplex_return(IT_0084);
}

