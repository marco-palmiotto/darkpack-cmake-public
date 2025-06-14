#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_m_se_R.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_m_se_R(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t M_eR = param->M_eR;
    const creal_t beta = param->beta;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = pow(M_eR, 2);
    const ccomplex_t IT_0001 = pow(m_e, 2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -2);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = pow(m_Z, 2);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = cpow(IT_0000 + IT_0001*IT_0003 + -IT_0001
      *IT_0003*IT_0005 + (-0.5)*IT_0006*IT_0008 + IT_0005*IT_0006*IT_0008, 0.5);
    return create_ccomplex_return(IT_0009);
}

