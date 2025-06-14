#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_m_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_m_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_A0 = param->m_A0;
    const ccomplex_t IT_0000 = pow(m_W, 2);
    const ccomplex_t IT_0001 = pow(m_A0, 2);
    const ccomplex_t IT_0002 = cpow(IT_0000 + IT_0001, 0.5);
    return create_ccomplex_return(IT_0002);
}

