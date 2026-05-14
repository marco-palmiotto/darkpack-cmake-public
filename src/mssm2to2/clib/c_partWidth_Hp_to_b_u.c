#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_Hp_to_b_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_Hp_to_b_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t theta_W = param->theta_W;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_Hp, -2);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_b*e_em
      *IT_0002*IT_0003*IT_0005*IT_0006*IT_0008*V_ub_mod;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0010, 2);
    const ccomplex_t IT_0012 = pow(m_u, 2);
    const ccomplex_t IT_0013 = pow(m_b, 2);
    const ccomplex_t IT_0014 = pow(m_Hp, 2);
    const ccomplex_t IT_0015 = cpow(-IT_0012 + IT_0013 + -IT_0014, 2);
    const ccomplex_t IT_0016 = -IT_0012 + (-0.25)*IT_0001*(IT_0012 + -IT_0013 
      + IT_0014)*(IT_0012 + -IT_0013 + -IT_0014) + 0.25*IT_0001*IT_0015;
    const ccomplex_t IT_0017 = 6*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0006, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0002*IT_0003*IT_0004*IT_0008*IT_0018*V_ub_mod;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = m_b*m_u;
    const ccomplex_t IT_0022 = cabs(-IT_0012 + 0.25*IT_0001*IT_0015);
    const ccomplex_t IT_0023 = cpow(IT_0022, 0.5);
    const ccomplex_t IT_0024 = ( ((m_Hp
) >= (m_b + m_u
)) ? (0.125*IT_0000
      *IT_0001*(IT_0011*IT_0017 + IT_0020*(IT_0017*IT_0020 + (-12)*IT_0010
      *IT_0021))*IT_0023
) : (0
));
    return create_ccomplex_return(IT_0024);
}

