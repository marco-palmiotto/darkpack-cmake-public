#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_H0_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_H0_to_h_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t alpha = param->alpha;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -2);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = cpow(IT_0002, -1);
    const ccomplex_t IT_0005 = IT_0003 + IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = cpow(IT_0007, 3);
    const ccomplex_t IT_0009 = cos(alpha);
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = cpow(IT_0009, 2);
    const ccomplex_t IT_0012 = cpow(IT_0007, 2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*m_W*e_em*IT_0005*(IT_0006
      *IT_0008 + (-5)*IT_0009*(IT_0006*IT_0007*IT_0009 + (-0.2)*IT_0010*(IT_0011
       + (-5)*IT_0012)));
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0014, 2);
    return create_ccomplex_return(-IT_0015);
}

