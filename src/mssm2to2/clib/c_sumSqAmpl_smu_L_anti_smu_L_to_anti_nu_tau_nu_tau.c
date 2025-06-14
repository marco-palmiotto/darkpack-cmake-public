#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_smu_L_to_anti_nu_tau_nu_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_smu_L_to_anti_nu_tau_nu_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0002;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0003;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0007 + IT_0009);
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = pow(m_smu_L, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = s_13*s_14;
    const ccomplex_t IT_0019 = s_34*IT_0013;
    return create_ccomplex_return(4*IT_0017*conj(IT_0017)*(IT_0018 + (-0.5)
      *IT_0019));
}

