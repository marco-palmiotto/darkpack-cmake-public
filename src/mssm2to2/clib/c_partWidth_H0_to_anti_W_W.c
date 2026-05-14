#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_H0_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_H0_to_anti_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t alpha = param->alpha;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_H0, -2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*m_W*e_em*(IT_0002*IT_0003 +
       IT_0004*IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = pow(m_W, -4);
    const ccomplex_t IT_0011 = pow(m_W, 2);
    const ccomplex_t IT_0012 = pow(m_H0, 2);
    const ccomplex_t IT_0013 = -IT_0011 + 0.5*IT_0012;
    const ccomplex_t IT_0014 = cpow(IT_0013, 2);
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = cabs(-IT_0011 + 0.25*IT_0012);
    const ccomplex_t IT_0017 = cpow(IT_0016, 0.5);
    const ccomplex_t IT_0018 = pow(m_W, -2);
    const ccomplex_t IT_0019 = pow(m_W, 4);
    const ccomplex_t IT_0020 = pow(m_H0, -4);
    const ccomplex_t IT_0021 = clog(IT_0001*IT_0011);
    const ccomplex_t IT_0022 = cpow(-1 + 4*IT_0001*IT_0011, (-0.5));
    const ccomplex_t IT_0023 = pow(m_W, -3);
    const ccomplex_t IT_0024 = pow(m_H0, 3);
    const ccomplex_t IT_0025 = cacos((-0.5)*(1 + (-3)*IT_0001*IT_0011)*IT_0023
      *IT_0024);
    const ccomplex_t IT_0026 = (1 + -IT_0001*IT_0011)*IT_0012*IT_0018*(1 + (
      -6.5)*IT_0001*IT_0011 + 23.5*IT_0019*IT_0020) + 1.5*(1 + (-6)*IT_0001
      *IT_0011 + 4*IT_0019*IT_0020)*IT_0021 + (-3)*(1 + (-8)*IT_0001*IT_0011 +
       20*IT_0019*IT_0020)*IT_0022*IT_0025;
    const ccomplex_t IT_0027 = ( ((m_H0
) >= (2*m_W
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0009*(1 + 0.5*IT_0015)*IT_0017
) : (0.25*IT_0000*IT_0001
      *IT_0009*(1 + 0.5*IT_0015)*IT_0017*IT_0026
));
    return create_ccomplex_return(IT_0027);
}

