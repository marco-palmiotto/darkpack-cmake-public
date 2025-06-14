#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_smu_R_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_smu_R_to_h_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0000
      *IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = pow(m_mu, 2);
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0005*IT_0008*IT_0009 + (-0.5)*m_W*(IT_0001*IT_0009 + IT_0003*IT_0010)*
      (IT_0005 + -IT_0004*IT_0012));
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = pow(m_h, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0015 + (0 + _Complex_I*1)
      *m_smu_L*Gamma_ml + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*-1)*e_em*(IT_0000*IT_0002
      *IT_0005*IT_0008*IT_0009 + -m_W*IT_0004*(IT_0001*IT_0009 + IT_0003*IT_0010
      )*IT_0012);
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_23 + IT_0015 + (0 + _Complex_I*1)
      *m_smu_R*Gamma_mr + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = pow(m_smu_R, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0015 + IT_0025 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0014*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = pow(m_smu_L, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0015 + IT_0029 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0021*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0000
      *IT_0002*IT_0005*IT_0009;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0004*IT_0012;
    const ccomplex_t IT_0036 = IT_0005 + IT_0035;
    const ccomplex_t IT_0037 = cpow(IT_0009, 3);
    const ccomplex_t IT_0038 = cpow(IT_0003, 2);
    const ccomplex_t IT_0039 = cpow(IT_0009, 2);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_W*e_em*IT_0036*(IT_0010
      *IT_0037 + (-5)*IT_0003*(IT_0003*IT_0009*IT_0010 + (-0.2)*IT_0001*(IT_0038
       + (-5)*IT_0039)));
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0034*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + -IT_0025 + -IT_0029 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*m_W*e_em*IT_0036*(IT_0001
      *IT_0037 + -IT_0003*(IT_0003*(IT_0001*IT_0009 + IT_0003*IT_0010) + 
      -IT_0010*IT_0039));
    const ccomplex_t IT_0047 = 1.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0007*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0025 + -IT_0029 + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0018 + IT_0024 + IT_0028 + IT_0032 + 
      -IT_0045 + -IT_0051;
    return create_ccomplex_return(IT_0052*conj(IT_0052));
}

