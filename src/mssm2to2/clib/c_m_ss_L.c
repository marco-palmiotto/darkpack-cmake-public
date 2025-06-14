#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_m_ss_L.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_m_ss_L(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_s = param->m_s;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t M_q1L = param->M_q1L;
    const creal_t M_q3L = param->M_q3L;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t IT_0000 = pow(M_q1L, 2);
    const ccomplex_t IT_0001 = pow(V_us, 2);
    const ccomplex_t IT_0002 = pow(M_q3L, 2);
    const ccomplex_t IT_0003 = pow(m_s, 2);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -2);
    const ccomplex_t IT_0006 = pow(m_W, 2);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = pow(m_Z, 2);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, 2);
    const ccomplex_t IT_0012 = V_cs*conj(V_cs)*IT_0000 + IT_0000*IT_0001 +
       V_ts*conj(V_ts)*IT_0002 + IT_0003*IT_0005 + -IT_0003*IT_0005*IT_0008 + (
      -0.5)*IT_0006*IT_0008 + (-0.0833333333333333)*IT_0008*IT_0009*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0012, 0.5);
    return create_ccomplex_return(IT_0013);
}

