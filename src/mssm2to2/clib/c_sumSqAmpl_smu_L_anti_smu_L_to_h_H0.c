#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_smu_L_to_h_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_smu_L_to_h_H0(
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
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(alpha);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = pow(m_W, -2);
    const ccomplex_t IT_0006 = pow(m_mu, 2);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -2);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0000*IT_0001*IT_0002*
      (IT_0004 + -(1 + -IT_0005*IT_0006*IT_0008)*IT_0010);
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = cpow(IT_0007, -1);
    const ccomplex_t IT_0015 = cpow(IT_0009, -1);
    const ccomplex_t IT_0016 = sin(beta);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0006
      *IT_0013*IT_0014*IT_0015 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0015)*(IT_0001
      *IT_0007 + -IT_0002*IT_0016));
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0004*IT_0009;
    const ccomplex_t IT_0020 = IT_0015 + IT_0019;
    const ccomplex_t IT_0021 = cpow(IT_0002, 3);
    const ccomplex_t IT_0022 = cpow(IT_0002, 2);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_W*e_em*IT_0020*(IT_0007
      *IT_0021 + 5*IT_0001*((-1)*IT_0001*(IT_0002*IT_0007 + 1./5*IT_0001*IT_0016
      ) + IT_0016*IT_0022));
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0018*IT_0024;
    const ccomplex_t IT_0026 = pow(m_smu_L, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + (-2)*IT_0026 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0006
      *IT_0013*IT_0014*IT_0015 + 0.5*m_W*(IT_0004*IT_0009 + -IT_0015)*(IT_0002
      *IT_0007 + IT_0001*IT_0016));
    const ccomplex_t IT_0031 = cpow(IT_0001, 2);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*m_W*e_em*IT_0020*(IT_0016
      *IT_0021 + (-5)*IT_0001*(IT_0001*IT_0002*IT_0016 + 1*IT_0007*(IT_0022 + 
      -1./5*IT_0031)));
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0026 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0002
      *IT_0013*IT_0014*IT_0015;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*m_mu*mu_h*IT_0001
      *IT_0013*IT_0014*IT_0015;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = pow(m_h, 2);
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0026 + IT_0043 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0018*IT_0030;
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0043 + (0 + _Complex_I*1)
      *m_smu_L*Gamma_ml + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0043 + (0 + _Complex_I*1)
      *m_smu_L*Gamma_ml + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0047*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_23 + IT_0026 + IT_0043 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0042*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = IT_0012 + -IT_0029 + -IT_0037 + IT_0046 +
       IT_0050 + IT_0053 + IT_0056;
    return create_ccomplex_return(IT_0057*conj(IT_0057));
}

