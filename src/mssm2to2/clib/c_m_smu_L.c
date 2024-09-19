#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_m_smu_L.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_m_smu_L(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t M_eL = param->M_eL;
    const creal_t beta = param->beta;
    const creal_t m_mu = param->m_mu;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = pow(m_W, 2);
    const ccomplex_t IT_0001 = pow(M_eL, 2);
    const ccomplex_t IT_0002 = pow(m_mu, 2);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = pow(m_Z, 2);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = cpow((-0.5)*IT_0000 + IT_0001 + IT_0002*IT_0004
       + -IT_0002*IT_0004*IT_0006 + (-0.25)*IT_0007*IT_0009 + 0.5*IT_0006
      *IT_0007*IT_0009, 0.5);
    return create_ccomplex_return(IT_0010);
}

