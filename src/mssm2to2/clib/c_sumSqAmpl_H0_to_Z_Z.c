#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_H0_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_H0_to_Z_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -4);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0000*IT_0001*IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = IT_0003*IT_0004*IT_0006*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0004, -3);
    const ccomplex_t IT_0010 = IT_0000*IT_0001*IT_0009;
    const ccomplex_t IT_0011 = IT_0006*IT_0007*IT_0009;
    const ccomplex_t IT_0012 = cpow(IT_0002, -2);
    const ccomplex_t IT_0013 = cpow(IT_0004, -1);
    const ccomplex_t IT_0014 = IT_0000*IT_0001*IT_0012*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0007*IT_0012*IT_0013;
    const ccomplex_t IT_0016 = cpow(IT_0004, -2);
    const ccomplex_t IT_0017 = cpow(IT_0012 + IT_0016, -1);
    const ccomplex_t IT_0018 = (IT_0005 + IT_0008 + IT_0010 + IT_0011 + 2
      *IT_0014 + 2*IT_0015)*IT_0017;
    const ccomplex_t IT_0019 = m_W*e_em;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = cpow(IT_0021, 2);
    const ccomplex_t IT_0023 = pow(m_Z, -4);
    const ccomplex_t IT_0024 = pow(s_23, 2);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    return create_ccomplex_return((-2)*IT_0022*(1 + 0.5*IT_0025));
}

