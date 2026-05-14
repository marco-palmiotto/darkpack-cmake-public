#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_smu_R_to_mu_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_smu_R_to_mu_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0002*IT_0004*IT_0006;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = pow(m_mu, 2);
    const ccomplex_t IT_0010 = pow(m_N_4, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = IT_0001*IT_0013;
    const ccomplex_t IT_0015 = pow(m_smu_R, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + IT_0009 + -IT_0010 + -IT_0015 
      + -reg_prop, -1);
    const ccomplex_t IT_0017 = -IT_0000;
    const ccomplex_t IT_0018 = IT_0008*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0014 + 2*IT_0020;
    const ccomplex_t IT_0022 = cos(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0011*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0001*IT_0027;
    const ccomplex_t IT_0029 = IT_0017*IT_0025;
    const ccomplex_t IT_0030 = IT_0016*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = -IT_0028 + (-2)*IT_0031;
    const ccomplex_t IT_0033 = (-2)*IT_0009;
    const ccomplex_t IT_0034 = IT_0000*IT_0013;
    const ccomplex_t IT_0035 = IT_0000*IT_0027;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = (-2)*s_13;
    const ccomplex_t IT_0038 = s_34*IT_0037;
    const ccomplex_t IT_0039 = -IT_0014;
    const ccomplex_t IT_0040 = -IT_0031;
    const ccomplex_t IT_0041 = m_mu*s_14*m_N_4;
    const ccomplex_t IT_0042 = s_13*s_14;
    const ccomplex_t IT_0043 = s_34*IT_0010;
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0042 + IT_0044;
    const ccomplex_t IT_0046 = s_14*IT_0037;
    const ccomplex_t IT_0047 = IT_0010*IT_0037;
    const ccomplex_t IT_0048 = pow(s_13, 2);
    const ccomplex_t IT_0049 = (-2)*IT_0048;
    const ccomplex_t IT_0050 = m_mu*m_N_4;
    const ccomplex_t IT_0051 = (-2)*s_14;
    const ccomplex_t IT_0052 = (-2)*s_34;
    const ccomplex_t IT_0053 = (-2)*IT_0010;
    const ccomplex_t IT_0054 = (conj(IT_0021)*IT_0032 + IT_0021*conj(IT_0032))
      *IT_0033 + IT_0037*(conj(IT_0021)*IT_0028 + IT_0021*conj(IT_0028) + conj
      (IT_0032)*IT_0039 + IT_0032*conj(IT_0039)) + (conj(IT_0028)*IT_0034 +
       IT_0028*conj(IT_0034) + conj(IT_0036)*IT_0039 + IT_0036*conj(IT_0039))
      *IT_0051 + (conj(IT_0032)*IT_0034 + IT_0032*conj(IT_0034) + conj(IT_0021)
      *IT_0036 + IT_0021*conj(IT_0036))*IT_0052 + (conj(IT_0028)*IT_0039 +
       IT_0028*conj(IT_0039))*IT_0053;
    const ccomplex_t IT_0055 = m_mu*s_34*m_N_4;
    const ccomplex_t IT_0056 = (-2)*conj(IT_0021);
    const ccomplex_t IT_0057 = s_14*IT_0009;
    const ccomplex_t IT_0058 = (-2)*conj(IT_0032);
    const ccomplex_t IT_0059 = s_14*s_34;
    const ccomplex_t IT_0060 = (-2)*IT_0032;
    const ccomplex_t IT_0061 = (-0.5)*IT_0057;
    const ccomplex_t IT_0062 = (-2)*IT_0021;
    const ccomplex_t IT_0063 = (-0.5)*IT_0055;
    const ccomplex_t IT_0064 = s_13*(IT_0021*conj(IT_0021) + IT_0032*conj
      (IT_0032))*IT_0033 + (conj(IT_0021)*IT_0034 + IT_0021*conj(IT_0034) + conj
      (IT_0032)*IT_0036 + IT_0032*conj(IT_0036))*IT_0038 + (-2)*conj(IT_0039)*
      (IT_0040*IT_0041 + 2*IT_0020*IT_0045 + (-0.5)*IT_0034*IT_0046 + (-0.5)
      *IT_0039*IT_0047 + (-0.5)*IT_0021*IT_0049) + IT_0039*(conj(IT_0034)
      *IT_0046 + conj(IT_0021)*IT_0049) + (-2)*conj(IT_0028)*(IT_0020*IT_0041 +
       2*IT_0040*IT_0045 + (-0.5)*IT_0036*IT_0046 + (-0.5)*IT_0028*IT_0047 + (
      -0.5)*IT_0032*IT_0049) + IT_0028*(conj(IT_0036)*IT_0046 + conj(IT_0032)
      *IT_0049) + IT_0050*IT_0054 + IT_0020*(IT_0056*IT_0057 + IT_0055*IT_0058) 
      + IT_0040*(IT_0055*IT_0056 + IT_0057*IT_0058) + (-2)*conj(IT_0020)*
      (IT_0028*IT_0041 + 2*IT_0039*IT_0045 + (-4)*IT_0020*IT_0059 + IT_0061
      *IT_0062 + IT_0060*IT_0063) + (-2)*conj(IT_0040)*(IT_0039*IT_0041 + 2
      *IT_0028*IT_0045 + (-4)*IT_0040*IT_0059 + IT_0060*IT_0061 + IT_0062
      *IT_0063);
    return create_ccomplex_return(IT_0064);
}

