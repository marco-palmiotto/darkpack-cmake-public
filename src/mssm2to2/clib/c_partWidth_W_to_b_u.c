#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_W_to_b_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_W_to_b_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_u = param->m_u;
    const creal_t e_em = param->e_em;
    const creal_t theta_W = param->theta_W;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0002
      *IT_0004*V_ub_mod;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = pow(m_u, 2);
    const ccomplex_t IT_0009 = pow(m_W, 2);
    const ccomplex_t IT_0010 = pow(m_b, 2);
    const ccomplex_t IT_0011 = cpow(IT_0008 + IT_0009 + -IT_0010, 2);
    const ccomplex_t IT_0012 = -IT_0008 + (-0.25)*IT_0001*(IT_0008 + IT_0009 +
       -IT_0010)*(IT_0008 + -IT_0009 + -IT_0010) + 0.25*IT_0001*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0008 + 0.5*IT_0009 + (-0.5)*IT_0010;
    const ccomplex_t IT_0014 = (-0.5)*IT_0008 + 0.5*IT_0009 + 0.5*IT_0010;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = IT_0009*IT_0012;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0015 + IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = cabs(-IT_0008 + 0.25*IT_0001*IT_0011);
    const ccomplex_t IT_0021 = cpow(IT_0020, 0.5);
    const ccomplex_t IT_0022 = ( ((m_W
) >= (m_b + m_u
)) ? (0.5*IT_0000
      *IT_0001*IT_0007*(IT_0012 + IT_0019)*IT_0021
) : (0
));
    return create_ccomplex_return(IT_0022);
}

