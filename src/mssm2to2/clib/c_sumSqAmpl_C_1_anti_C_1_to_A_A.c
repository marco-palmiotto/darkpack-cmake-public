#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_A_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_A_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_C_1
      *Gamma_c1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*IT_0002;
    const ccomplex_t IT_0004 = IT_0000*IT_0003;
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0004;
    const ccomplex_t IT_0007 = pow(m_C_1, 2);
    const ccomplex_t IT_0008 = 16*s_14 + 8*s_24;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_1
      *Gamma_c1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = s_12*s_14;
    const ccomplex_t IT_0014 = s_24*IT_0007;
    const ccomplex_t IT_0015 = 8*IT_0014;
    const ccomplex_t IT_0016 = s_14*IT_0007;
    const ccomplex_t IT_0017 = (-16)*IT_0016;
    const ccomplex_t IT_0018 = (-16)*IT_0013 + IT_0015 + IT_0017;
    const ccomplex_t IT_0019 = 2*IT_0007;
    const ccomplex_t IT_0020 = s_12 + IT_0019;
    const ccomplex_t IT_0021 = (-8)*IT_0007;
    const ccomplex_t IT_0022 = 2*s_24;
    const ccomplex_t IT_0023 = s_14 + IT_0022;
    const ccomplex_t IT_0024 = IT_0007*IT_0023;
    const ccomplex_t IT_0025 = s_12*s_24;
    const ccomplex_t IT_0026 = (-2)*IT_0025;
    const ccomplex_t IT_0027 = (-2)*IT_0014 + IT_0016 + IT_0026;
    const ccomplex_t IT_0028 = s_14*s_24;
    const ccomplex_t IT_0029 = (-8)*IT_0011;
    const ccomplex_t IT_0030 = -conj(IT_0005)*(IT_0006*IT_0009 + -IT_0012
      *IT_0018 + 8*IT_0020*(s_12*IT_0011 + (-0.125)*IT_0005*IT_0021)) + (-8)
      *conj(IT_0011)*(IT_0020*(s_12*IT_0005 + (-0.125)*IT_0011*IT_0021) +
       IT_0012*IT_0024 + -IT_0006*IT_0027) + conj(IT_0012)*(IT_0005*IT_0018 + 32
      *IT_0012*IT_0028 + IT_0024*IT_0029) + -conj(IT_0006)*(IT_0005*IT_0009 + (
      -32)*IT_0006*IT_0028 + IT_0027*IT_0029);
    return create_ccomplex_return(IT_0030);
}

