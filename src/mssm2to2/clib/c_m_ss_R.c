#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_m_ss_R.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_m_ss_R(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t m_s = param->m_s;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t M_qbR = param->M_qbR;
    const creal_t M_qdR = param->M_qdR;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t IT_0000 = pow(M_qbR, 2);
    const ccomplex_t IT_0001 = pow(M_qdR, 2);
    const ccomplex_t IT_0002 = pow(V_us, 2);
    const ccomplex_t IT_0003 = pow(m_s, 2);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -2);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = pow(m_Z, 2);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, 2);
    const ccomplex_t IT_0011 = V_ts*conj(V_ts)*IT_0000 + V_cs*conj(V_cs)
      *IT_0001 + IT_0001*IT_0002 + IT_0003*IT_0005 + -IT_0003*IT_0005*IT_0007 + 
      (-0.166666666666667)*IT_0008*IT_0010 + 0.333333333333333*IT_0007*IT_0008
      *IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0011, 0.5);
    return create_ccomplex_return(IT_0012);
}

