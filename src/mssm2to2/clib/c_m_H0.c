#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_m_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_m_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t m_A0 = param->m_A0;
    const ccomplex_t IT_0000 = pow(m_Z, 2);
    const ccomplex_t IT_0001 = pow(m_A0, 2);
    const ccomplex_t IT_0002 = cpow(-IT_0000 + IT_0001, 2);
    const ccomplex_t IT_0003 = sin(2*beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = cpow(IT_0002 + 4*IT_0000*IT_0001*IT_0004, 0.5);
    const ccomplex_t IT_0006 = cpow(0.5*IT_0000 + 0.5*IT_0001 + 0.5*IT_0005,
       0.5);
    return create_ccomplex_return(IT_0006);
}

