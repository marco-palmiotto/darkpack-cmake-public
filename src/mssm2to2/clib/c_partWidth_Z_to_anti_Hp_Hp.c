#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_partWidth_Z_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_partWidth_Z_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0004 + 
      -IT_0002*IT_0005);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = pow(m_Hp, 2);
    const ccomplex_t IT_0009 = pow(m_Z, 2);
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0010, 2);
    const ccomplex_t IT_0012 = IT_0001*IT_0011;
    const ccomplex_t IT_0013 = cabs(-IT_0008 + 0.25*IT_0009);
    const ccomplex_t IT_0014 = cpow(IT_0013, 0.5);
    const ccomplex_t IT_0015 = ( ((m_Z
) >= (2*m_Hp
)) ? (0.0416666666666667
      *IT_0000*IT_0001*IT_0007*(IT_0008 + -IT_0012)*IT_0014
) : (0
));
    return create_ccomplex_return(IT_0015);
}

