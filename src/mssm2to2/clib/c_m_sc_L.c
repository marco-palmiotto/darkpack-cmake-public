#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_m_sc_L.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_m_sc_L(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t beta = param->beta;
    const creal_t M_q1L = param->M_q1L;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = pow(m_W, 2);
    const ccomplex_t IT_0001 = pow(m_c, 2);
    const ccomplex_t IT_0002 = pow(M_q1L, 2);
    const ccomplex_t IT_0003 = pow(m_Z, 2);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = cpow((-0.5)*IT_0000 + IT_0001 + IT_0002 +
       0.0833333333333333*IT_0003*IT_0005 + (-0.166666666666667)*IT_0003*IT_0005
      *IT_0007, 0.5);
    return create_ccomplex_return(IT_0008);
}

