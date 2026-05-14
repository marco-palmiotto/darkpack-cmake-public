#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_smu_R_to_H0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_smu_R_to_H0_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = e_em*m_mu*mu_h*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = IT_0002 + IT_0007;
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, 3);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = cos(alpha);
    const ccomplex_t IT_0015 = cpow(IT_0009, 2);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*m_W*e_em*IT_0008*(IT_0010
      *IT_0011 + -IT_0013*(IT_0013*(IT_0009*IT_0011 + -IT_0012*IT_0014) + 
      -IT_0014*IT_0015));
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = IT_0004*IT_0017;
    const ccomplex_t IT_0019 = pow(m_smu_L, 2);
    const ccomplex_t IT_0020 = pow(m_smu_R, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + -IT_0019 + -IT_0020 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = pow(m_mu, 2);
    const ccomplex_t IT_0025 = cpow(IT_0012, -1);
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*e_em*(0.5*m_W*(IT_0002 + 
      -IT_0001*IT_0006)*(IT_0009*IT_0011 + -IT_0012*IT_0014) + IT_0000*IT_0002
      *IT_0014*IT_0024*IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0004*IT_0027;
    const ccomplex_t IT_0029 = pow(m_H0, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0029 + (0 + _Complex_I*1)
      *m_smu_L*Gamma_ml + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*-1)*e_em*(m_W*IT_0001*IT_0006*
      (IT_0009*IT_0011 + -IT_0012*IT_0014) + IT_0000*IT_0002*IT_0014*IT_0024
      *IT_0025);
    const ccomplex_t IT_0034 = IT_0004*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0029 + (0 + _Complex_I*1)
      *m_smu_R*Gamma_mr + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*-1)*m_W*e_em*IT_0008*IT_0009
      *IT_0012*(IT_0009*IT_0011 + -IT_0012*IT_0014);
    const ccomplex_t IT_0039 = e_em*m_mu*mu_h*IT_0000*IT_0002*IT_0009*IT_0025;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = pow(m_Z, 2);
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + -IT_0019 + -IT_0020 + IT_0042 
      + -reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = -IT_0023 + IT_0032 + IT_0037 + -IT_0045;
    return create_ccomplex_return(IT_0046*conj(IT_0046));
}

