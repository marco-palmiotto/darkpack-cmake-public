#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_m_chi_dm_1.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_m_chi_dm_1(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_cdm = param->m_cdm;
    const ccomplex_t IT_0000 = (-0.5)*m_cdm;
    return create_ccomplex_return(IT_0000);
}

