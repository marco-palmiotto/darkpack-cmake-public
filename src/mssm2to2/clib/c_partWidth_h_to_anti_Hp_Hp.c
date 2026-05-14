#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_h_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_h_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t alpha = param->alpha;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_h, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -2);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0004, -1);
    const ccomplex_t IT_0007 = IT_0005 + -IT_0006;
    const ccomplex_t IT_0008 = cos(alpha);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, 3);
    const ccomplex_t IT_0011 = cos(beta);
    const ccomplex_t IT_0012 = IT_0005 + 3*IT_0006;
    const ccomplex_t IT_0013 = cpow(IT_0009, 2);
    const ccomplex_t IT_0014 = sin(alpha);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_W*e_em*(IT_0007*IT_0008
      *IT_0010 + IT_0011*(IT_0012*IT_0013*IT_0014 + -IT_0011*(IT_0008*IT_0009
      *IT_0012 + IT_0007*IT_0011*IT_0014)));
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = cpow(IT_0016, 2);
    const ccomplex_t IT_0018 = pow(m_h, 2);
    const ccomplex_t IT_0019 = pow(m_Hp, 2);
    const ccomplex_t IT_0020 = cabs(0.25*IT_0018 + -IT_0019);
    const ccomplex_t IT_0021 = cpow(IT_0020, 0.5);
    const ccomplex_t IT_0022 = ( ((m_h
) >= (2*m_Hp
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0017*IT_0021
) : (0
));
    return create_ccomplex_return(IT_0022);
}

