#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_anti_ss_R_to_Hp_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_anti_ss_R_to_Hp_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 +
       IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = pow(m_s, 2);
    const ccomplex_t IT_0012 = cpow(IT_0000, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(m_W*(IT_0001*IT_0002 +
       IT_0000*IT_0003)*IT_0004*IT_0009 + (-3)*IT_0003*IT_0005*IT_0010*IT_0011
      *IT_0012);
    const ccomplex_t IT_0014 = (-0.333333333333333)*IT_0013;
    const ccomplex_t IT_0015 = pow(m_ss_R, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0007*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003)*IT_0005;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(m_W*(IT_0000*IT_0001 +
       -IT_0002*IT_0003)*IT_0004*IT_0009 + (-3)*IT_0001*IT_0005*IT_0010*IT_0011
      *IT_0012);
    const ccomplex_t IT_0022 = 0.333333333333333*IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + (-2)*IT_0015 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = IT_0019 + IT_0026;
    const ccomplex_t IT_0028 = pow(m_Hp, 2);
    const ccomplex_t IT_0029 = pow(m_W, -2);
    const ccomplex_t IT_0030 = pow(s_34, 2);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    return create_ccomplex_return((-3)*IT_0027*conj(IT_0027)*(IT_0028 + 
      -IT_0031));
}

