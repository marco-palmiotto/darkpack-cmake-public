#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_Hp_to_c_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_Hp_to_c_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_Hp, -2);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0002*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0003, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0002*IT_0004*IT_0007*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = m_c*m_s;
    const ccomplex_t IT_0014 = pow(m_s, 2);
    const ccomplex_t IT_0015 = pow(m_c, 2);
    const ccomplex_t IT_0016 = pow(m_Hp, 2);
    const ccomplex_t IT_0017 = cpow(-IT_0014 + IT_0015 + -IT_0016, 2);
    const ccomplex_t IT_0018 = -IT_0014 + (-0.25)*IT_0001*(IT_0014 + -IT_0015 
      + IT_0016)*(IT_0014 + -IT_0015 + -IT_0016) + 0.25*IT_0001*IT_0017;
    const ccomplex_t IT_0019 = 6*IT_0018;
    const ccomplex_t IT_0020 = cabs(-IT_0014 + 0.25*IT_0001*IT_0017);
    const ccomplex_t IT_0021 = cpow(IT_0020, 0.5);
    const ccomplex_t IT_0022 = ( ((m_Hp
) >= (m_c + m_s
)) ? ((-0.75)*IT_0000
      *IT_0001*(conj(IT_0009)*(IT_0012*IT_0013 + (-0.166666666666667)*IT_0009
      *IT_0019) + conj(IT_0012)*(IT_0009*IT_0013 + (-0.166666666666667)*IT_0012
      *IT_0019))*IT_0021
) : (0
));
    return create_ccomplex_return(IT_0022);
}

