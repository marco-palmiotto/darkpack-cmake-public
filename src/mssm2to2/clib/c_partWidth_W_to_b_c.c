#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_W_to_b_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_W_to_b_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t V_cb = param->V_cb;
    const creal_t e_em = param->e_em;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0003;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = pow(m_c, 2);
    const ccomplex_t IT_0008 = pow(m_W, 2);
    const ccomplex_t IT_0009 = pow(m_b, 2);
    const ccomplex_t IT_0010 = cpow(IT_0007 + IT_0008 + -IT_0009, 2);
    const ccomplex_t IT_0011 = -IT_0007 + (-0.25)*IT_0001*(IT_0007 + IT_0008 +
       -IT_0009)*(IT_0007 + -IT_0008 + -IT_0009) + 0.25*IT_0001*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0007 + 0.5*IT_0008 + (-0.5)*IT_0009;
    const ccomplex_t IT_0013 = (-0.5)*IT_0007 + 0.5*IT_0008 + 0.5*IT_0009;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = IT_0008*IT_0011;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0014 + IT_0016;
    const ccomplex_t IT_0018 = IT_0001*IT_0017;
    const ccomplex_t IT_0019 = cabs(-IT_0007 + 0.25*IT_0001*IT_0010);
    const ccomplex_t IT_0020 = cpow(IT_0019, 0.5);
    const ccomplex_t IT_0021 = ( ((m_W
) >= (m_b + m_c
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0006*(IT_0011 + IT_0018)*IT_0020
) : (0
));
    return create_ccomplex_return(IT_0021);
}

