#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_H0_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_H0_to_Z_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t alpha = param->alpha;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_H0, -2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -4);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = IT_0002*IT_0003*IT_0005*IT_0006;
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = IT_0005*IT_0006*IT_0008*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0006, -3);
    const ccomplex_t IT_0012 = IT_0002*IT_0003*IT_0011;
    const ccomplex_t IT_0013 = IT_0008*IT_0009*IT_0011;
    const ccomplex_t IT_0014 = cpow(IT_0004, -2);
    const ccomplex_t IT_0015 = cpow(IT_0006, -1);
    const ccomplex_t IT_0016 = IT_0002*IT_0003*IT_0014*IT_0015;
    const ccomplex_t IT_0017 = IT_0008*IT_0009*IT_0014*IT_0015;
    const ccomplex_t IT_0018 = cpow(IT_0006, -2);
    const ccomplex_t IT_0019 = cpow(IT_0014 + IT_0018, -1);
    const ccomplex_t IT_0020 = (IT_0007 + IT_0010 + IT_0012 + IT_0013 + 2
      *IT_0016 + 2*IT_0017)*IT_0019;
    const ccomplex_t IT_0021 = m_W*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0023, 2);
    const ccomplex_t IT_0025 = pow(m_Z, -4);
    const ccomplex_t IT_0026 = pow(m_Z, 2);
    const ccomplex_t IT_0027 = pow(m_H0, 2);
    const ccomplex_t IT_0028 = -IT_0026 + 0.5*IT_0027;
    const ccomplex_t IT_0029 = cpow(IT_0028, 2);
    const ccomplex_t IT_0030 = IT_0025*IT_0029;
    const ccomplex_t IT_0031 = cabs(-IT_0026 + 0.25*IT_0027);
    const ccomplex_t IT_0032 = cpow(IT_0031, 0.5);
    const ccomplex_t IT_0033 = pow(m_Z, -2);
    const ccomplex_t IT_0034 = pow(m_Z, 4);
    const ccomplex_t IT_0035 = pow(m_H0, -4);
    const ccomplex_t IT_0036 = clog(IT_0001*IT_0026);
    const ccomplex_t IT_0037 = cpow(-1 + 4*IT_0001*IT_0026, (-0.5));
    const ccomplex_t IT_0038 = pow(m_Z, -3);
    const ccomplex_t IT_0039 = pow(m_H0, 3);
    const ccomplex_t IT_0040 = cacos((-0.5)*(1 + (-3)*IT_0001*IT_0026)*IT_0038
      *IT_0039);
    const ccomplex_t IT_0041 = (1 + -IT_0001*IT_0026)*IT_0027*IT_0033*(1 + (
      -6.5)*IT_0001*IT_0026 + 23.5*IT_0034*IT_0035) + 1.5*(1 + (-6)*IT_0001
      *IT_0026 + 4*IT_0034*IT_0035)*IT_0036 + (-3)*(1 + (-8)*IT_0001*IT_0026 +
       20*IT_0034*IT_0035)*IT_0037*IT_0040;
    const ccomplex_t IT_0042 = ( ((m_H0
) >= (2*m_Z
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0024*(1 + 0.5*IT_0030)*IT_0032
) : (0.125*IT_0000*IT_0001
      *IT_0024*(1 + 0.5*IT_0030)*IT_0032*IT_0041
));
    return create_ccomplex_return(IT_0042);
}

