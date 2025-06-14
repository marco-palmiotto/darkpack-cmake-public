#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_anti_se_R_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_anti_se_R_to_Z_Z(
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
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = m_W*e_em;
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -4);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = IT_0001*IT_0003*IT_0004*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = IT_0003*IT_0005*IT_0007*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0005, -3);
    const ccomplex_t IT_0011 = IT_0001*IT_0004*IT_0010;
    const ccomplex_t IT_0012 = IT_0007*IT_0008*IT_0010;
    const ccomplex_t IT_0013 = cpow(IT_0002, -2);
    const ccomplex_t IT_0014 = cpow(IT_0005, -1);
    const ccomplex_t IT_0015 = IT_0001*IT_0004*IT_0013*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0008*IT_0013*IT_0014;
    const ccomplex_t IT_0017 = cpow(IT_0005, -2);
    const ccomplex_t IT_0018 = cpow(IT_0013 + IT_0017, -1);
    const ccomplex_t IT_0019 = (IT_0006 + -IT_0009 + IT_0011 + -IT_0012 + 2
      *IT_0015 + (-2)*IT_0016)*IT_0018;
    const ccomplex_t IT_0020 = IT_0000*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = pow(m_W, -1);
    const ccomplex_t IT_0023 = cpow(IT_0007, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0001
      *IT_0014*IT_0022*IT_0023;
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0021*IT_0025;
    const ccomplex_t IT_0027 = pow(m_se_L, 2);
    const ccomplex_t IT_0028 = pow(m_se_R, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0027 + -IT_0028 + -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0001*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0033 = IT_0003*IT_0004*IT_0005*IT_0008;
    const ccomplex_t IT_0034 = IT_0001*IT_0007*IT_0010;
    const ccomplex_t IT_0035 = IT_0004*IT_0008*IT_0010;
    const ccomplex_t IT_0036 = IT_0001*IT_0007*IT_0013*IT_0014;
    const ccomplex_t IT_0037 = IT_0004*IT_0008*IT_0013*IT_0014;
    const ccomplex_t IT_0038 = IT_0018*(IT_0032 + IT_0033 + IT_0034 + IT_0035 
      + 2*IT_0036 + 2*IT_0037);
    const ccomplex_t IT_0039 = IT_0000*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0008
      *IT_0014*IT_0022*IT_0023;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_12 + -IT_0027 + -IT_0028 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = -IT_0031 + -IT_0046;
    const ccomplex_t IT_0048 = pow(m_Z, -4);
    const ccomplex_t IT_0049 = pow(s_34, 2);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    return create_ccomplex_return(2*IT_0047*conj(IT_0047)*(1 + 0.5*IT_0050));
}

