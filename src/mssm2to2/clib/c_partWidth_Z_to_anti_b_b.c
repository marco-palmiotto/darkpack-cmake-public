#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_Z_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_Z_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t e_em = param->e_em;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0004;
    const ccomplex_t IT_0006 = 0.333333333333333*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = pow(m_Z, 2);
    const ccomplex_t IT_0009 = pow(m_b, 2);
    const ccomplex_t IT_0010 = 0.5*IT_0008 + -IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0008;
    const ccomplex_t IT_0012 = cpow(IT_0011, 2);
    const ccomplex_t IT_0013 = IT_0008*IT_0010;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0012 + IT_0014;
    const ccomplex_t IT_0016 = IT_0001*IT_0015;
    const ccomplex_t IT_0017 = IT_0010 + IT_0016;
    const ccomplex_t IT_0018 = IT_0003*IT_0004;
    const ccomplex_t IT_0019 = e_em*IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0004, -1);
    const ccomplex_t IT_0021 = IT_0002*IT_0020;
    const ccomplex_t IT_0022 = e_em*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0019 + 3*IT_0022);
    const ccomplex_t IT_0024 = (-0.166666666666667)*IT_0023;
    const ccomplex_t IT_0025 = cabs(0.25*IT_0008 + -IT_0009);
    const ccomplex_t IT_0026 = cpow(IT_0025, 0.5);
    const ccomplex_t IT_0027 = ( ((m_Z
) >= (2*m_b
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0007*IT_0017 + (-0.25)*IT_0024*((-12)*IT_0006*IT_0009 + (-4)
      *IT_0017*IT_0024))*IT_0026
) : (0
));
    return create_ccomplex_return(IT_0027);
}

