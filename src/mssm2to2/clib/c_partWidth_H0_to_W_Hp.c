#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_H0_to_W_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_H0_to_W_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t alpha = param->alpha;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_H0, -2);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0003 + 
      -IT_0004*IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = pow(m_H0, 2);
    const ccomplex_t IT_0011 = pow(m_W, -2);
    const ccomplex_t IT_0012 = pow(m_W, 2);
    const ccomplex_t IT_0013 = pow(m_Hp, 2);
    const ccomplex_t IT_0014 = 0.5*IT_0010 + 0.5*IT_0012 + (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0014, 2);
    const ccomplex_t IT_0016 = IT_0011*IT_0015;
    const ccomplex_t IT_0017 = cpow(IT_0010 + IT_0012 + -IT_0013, 2);
    const ccomplex_t IT_0018 = cabs(-IT_0012 + 0.25*IT_0001*IT_0017);
    const ccomplex_t IT_0019 = cpow(IT_0018, 0.5);
    const ccomplex_t IT_0020 = ( ((m_H0
) >= (m_W + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*IT_0009*(IT_0010 + -IT_0016)*IT_0019
) : (0
));
    return create_ccomplex_return(IT_0020);
}

