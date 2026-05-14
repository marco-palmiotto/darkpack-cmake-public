#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_t_to_s_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_t_to_s_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_t, -2);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0002*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0003, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0002*IT_0004*IT_0007*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = m_s*m_t;
    const ccomplex_t IT_0014 = pow(m_s, 2);
    const ccomplex_t IT_0015 = pow(m_t, 2);
    const ccomplex_t IT_0016 = pow(m_Hp, 2);
    const ccomplex_t IT_0017 = 0.5*IT_0014 + 0.5*IT_0015 + (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0014 + -IT_0015 + -IT_0016, 2);
    const ccomplex_t IT_0019 = cabs(-IT_0016 + 0.25*IT_0001*IT_0018);
    const ccomplex_t IT_0020 = cpow(IT_0019, 0.5);
    const ccomplex_t IT_0021 = ( ((m_t
) >= (m_s + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*((conj(IT_0009)*IT_0012 + IT_0009*conj(IT_0012))*IT_0013 + 
      (IT_0009*conj(IT_0009) + IT_0012*conj(IT_0012))*IT_0017)*IT_0020
) : (0
));
    return create_ccomplex_return(IT_0021);
}

