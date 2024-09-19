#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_su_R_to_A0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_su_R_to_A0_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t alpha = param->alpha;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0000
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = IT_0005 + IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0002, 3);
    const ccomplex_t IT_0013 = sin(alpha);
    const ccomplex_t IT_0014 = cos(beta);
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = cpow(IT_0002, 2);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*m_W*e_em*IT_0011*(IT_0012
      *IT_0013 + -IT_0015*((IT_0002*IT_0013 + -IT_0001*IT_0014)*IT_0015 + 
      -IT_0001*IT_0016));
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = IT_0007*IT_0018;
    const ccomplex_t IT_0020 = pow(m_su_L, 2);
    const ccomplex_t IT_0021 = pow(m_su_R, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + -IT_0020 + -IT_0021 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0000
      *IT_0003*IT_0005*IT_0013;
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*m_W*e_em*IT_0011*(IT_0001
      *IT_0012 + -IT_0014*(IT_0014*(IT_0001*IT_0002 + IT_0013*IT_0014) + 
      -IT_0013*IT_0016));
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0020 + -IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = -IT_0024 + -IT_0032;
    return create_ccomplex_return(3*IT_0033*conj(IT_0033));
}

