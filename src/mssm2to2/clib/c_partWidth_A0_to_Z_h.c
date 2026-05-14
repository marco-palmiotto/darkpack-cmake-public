#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_A0_to_Z_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_A0_to_Z_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t alpha = param->alpha;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_A0, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = e_em*IT_0003*(IT_0004*IT_0005 + IT_0006*IT_0007
      )*IT_0009;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0012, 2);
    const ccomplex_t IT_0014 = pow(m_A0, 2);
    const ccomplex_t IT_0015 = pow(m_Z, -2);
    const ccomplex_t IT_0016 = pow(m_Z, 2);
    const ccomplex_t IT_0017 = pow(m_h, 2);
    const ccomplex_t IT_0018 = 0.5*IT_0014 + 0.5*IT_0016 + (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = cpow(IT_0018, 2);
    const ccomplex_t IT_0020 = IT_0015*IT_0019;
    const ccomplex_t IT_0021 = cpow(IT_0014 + IT_0016 + -IT_0017, 2);
    const ccomplex_t IT_0022 = cabs(-IT_0016 + 0.25*IT_0001*IT_0021);
    const ccomplex_t IT_0023 = cpow(IT_0022, 0.5);
    const ccomplex_t IT_0024 = ( ((m_A0
) >= (m_Z + m_h
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0013*(IT_0014 + -IT_0020)*IT_0023
) : (0
));
    return create_ccomplex_return(IT_0024);
}

