#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_h_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_h_to_anti_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t e_em = param->e_em;
    const creal_t s_23 = param->s_23;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = pow(m_W, -4);
    const ccomplex_t IT_0001 = pow(s_23, 2);
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*m_W*e_em*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    return create_ccomplex_return((-2)*(1 + 0.5*IT_0002)*IT_0006);
}

