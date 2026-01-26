#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_anti_chi_chi_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_anti_chi_chi_to_anti_tau_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_l = param->g_l;
    const creal_t m_W = param->m_W;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_34 = param->s_34;
    const creal_t g_chi = param->g_chi;
    const creal_t m_chi = param->m_chi;
    const creal_t m_phi = param->m_phi;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_phi = param->Gamma_phi;
    const ccomplex_t IT_0000 = s_12*s_34;
    const ccomplex_t IT_0001 = pow(m_chi, 2);
    const ccomplex_t IT_0002 = s_34*IT_0001;
    const ccomplex_t IT_0003 = pow(m_tau, 2);
    const ccomplex_t IT_0004 = s_12*IT_0003;
    const ccomplex_t IT_0005 = IT_0001*IT_0003;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*g_chi;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*g_l*e_em*m_tau*IT_0008
      *IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = cpow(s_12 + IT_0001 + 0.5*reg_prop + (-0.5)
      *m_phi*(m_phi + (0 + _Complex_I*-1)*Gamma_phi), -1);
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    return create_ccomplex_return((16*IT_0000 + (-16)*IT_0002 + (-16)*IT_0004 
      + 16*IT_0005)*IT_0017*conj(IT_0017));
}

