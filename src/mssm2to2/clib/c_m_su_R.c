#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_m_su_R.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_m_su_R(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t M_quR = param->M_quR;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(M_quR, 2);
    const ccomplex_t IT_0002 = pow(m_Z, 2);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = cpow(IT_0000 + IT_0001 + 0.333333333333333
      *IT_0002*IT_0004 + (-0.666666666666667)*IT_0002*IT_0004*IT_0006, 0.5);
    return create_ccomplex_return(IT_0007);
}

