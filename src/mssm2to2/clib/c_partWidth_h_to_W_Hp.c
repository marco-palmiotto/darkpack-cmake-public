#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_h_to_W_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_h_to_W_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t alpha = param->alpha;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_h, -2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0003 +
       IT_0004*IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = -IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0009, 2);
    const ccomplex_t IT_0011 = pow(m_h, 2);
    const ccomplex_t IT_0012 = pow(m_W, -2);
    const ccomplex_t IT_0013 = pow(m_W, 2);
    const ccomplex_t IT_0014 = pow(m_Hp, 2);
    const ccomplex_t IT_0015 = 0.5*IT_0011 + 0.5*IT_0013 + (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = cpow(IT_0015, 2);
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0011 + IT_0013 + -IT_0014, 2);
    const ccomplex_t IT_0019 = cabs(-IT_0013 + 0.25*IT_0001*IT_0018);
    const ccomplex_t IT_0020 = cpow(IT_0019, 0.5);
    const ccomplex_t IT_0021 = ( ((m_h
) >= (m_W + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*IT_0010*(IT_0011 + -IT_0017)*IT_0020
) : (0
));
    return create_ccomplex_return(IT_0021);
}

