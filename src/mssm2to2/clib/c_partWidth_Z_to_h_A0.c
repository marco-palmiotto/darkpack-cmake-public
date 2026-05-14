#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_Z_to_h_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_Z_to_h_A0(
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
    const ccomplex_t IT_0001 = pow(m_Z, -2);
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
    const ccomplex_t IT_0012 = cpow(IT_0011, 2);
    const ccomplex_t IT_0013 = pow(m_A0, 2);
    const ccomplex_t IT_0014 = pow(m_Z, 2);
    const ccomplex_t IT_0015 = pow(m_h, 2);
    const ccomplex_t IT_0016 = 0.5*IT_0013 + 0.5*IT_0014 + (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = cpow(IT_0016, 2);
    const ccomplex_t IT_0018 = IT_0001*IT_0017;
    const ccomplex_t IT_0019 = cpow(IT_0013 + IT_0014 + -IT_0015, 2);
    const ccomplex_t IT_0020 = cabs(-IT_0013 + 0.25*IT_0001*IT_0019);
    const ccomplex_t IT_0021 = cpow(IT_0020, 0.5);
    const ccomplex_t IT_0022 = ( ((m_Z
) >= (m_h + m_A0
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0012*(IT_0013 + -IT_0018)*IT_0021

      ) : (0
));
    return create_ccomplex_return(IT_0022);
}

