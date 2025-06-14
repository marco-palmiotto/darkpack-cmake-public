#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_sd_R_to_anti_mu_nu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_sd_R_to_anti_mu_nu_mu(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cpow(IT_0003, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_d*m_u*V_ud
      *e_em*IT_0000*IT_0002*IT_0005*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = pow(m_sd_R, 2);
    const ccomplex_t IT_0013 = pow(m_su_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    return create_ccomplex_return(6*s_34*IT_0017*conj(IT_0017));
}

